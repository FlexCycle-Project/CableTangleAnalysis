/*
 * In-process bridge to the bundled KymoKnot library.
 *
 * This mirrors the pipeline in KymoKnot-clone/KymoKnot_linear.c without the
 * argv/file/process layer:
 *   build a linear KNTarc from memory -> qhull hybrid closure -> random
 *   rotation -> rectification -> knot identification -> TD/BU localization.
 *
 * KymoKnot and qhull keep process-global state, so analyses are serialized.
 * The intended execution model is one seed per process.
 */

#include "KymoKnotNativeBridge.h"

#include "KNT_arc.h"
#include "KNT_closures.h"
#include "KNT_defaults.h"
#include "KNT_identify.h"
#include "KNT_lib.h"
#include "KNT_qhull.h"
#include "KNT_simplify.h"
#include "lib/my_memory.h"
#include "lib/prng.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static KNTID_wspace *g_workspace = NULL;
static int g_workspace_alex_dim = 0;
static int g_workspace_buffer_dim = 0;
static unsigned int g_seed_counter = 0;

static KNTID_wspace *ensure_workspace(int n)
{
    const int required_alex = (2 * n > 2000) ? 2 * n : 2000;
    const int required_buffer = (n > 1000) ? n : 1000;

    if (g_workspace != NULL &&
        required_alex <= g_workspace_alex_dim &&
        required_buffer <= g_workspace_buffer_dim)
    {
        return g_workspace;
    }

    if (g_workspace != NULL)
    {
        KNTID_free_wspace(g_workspace);
        g_workspace = NULL;
    }

    KNTID_ws_var init;
    init.alex_dim = required_alex;
    init.buffer_dim = required_buffer;
    g_workspace = KNTID_wspace_init(&init);
    if (g_workspace == NULL)
        return NULL;

    g_workspace_alex_dim = required_alex;
    g_workspace_buffer_dim = required_buffer;
    return g_workspace;
}

static KNTarc *build_linear_arc(const double *xyz, int n)
{
    KNTarc *arc = KNTinit_ptr();
    if (arc == NULL)
        return NULL;

    KNTarc_SET_value(arc, len, n, int);
    KNTarc_SET_value(arc, start, 0, int);
    KNTarc_SET_value(arc, end, n - 1, int);
    KNTarc_SET_value(arc, arc_type, ARC_ID_LIN, char);
    strncpy(arc->closure, NONE, strlen(NONE));
    strncpy(arc->simplification, SIMPNONE, strlen(SIMPNONE));

    arc->coord = d2t(n, 3);
    arc->index = i1t(n);
    if (arc->coord == NULL || arc->index == NULL)
        return arc;

    for (int i = 0; i < n; ++i)
    {
        arc->coord[i][0] = xyz[3 * i + 0];
        arc->coord[i][1] = xyz[3 * i + 1];
        arc->coord[i][2] = xyz[3 * i + 2];
        arc->index[i] = i;
    }
    return arc;
}

static int get_idx_rect_chain_local(KNTarc *knt_rect, int idx_orig)
{
    for (int i = 0; i < knt_rect->len; ++i)
    {
        if (idx_orig == knt_rect->index[i])
            return i;
    }
    return -1;
}

static int append_segment(
    KymoKnotNativeResult *out,
    int *capacity,
    const char *knotType,
    int start,
    int end,
    int length)
{
    if (out->count >= *capacity)
    {
        int next = (*capacity == 0) ? 4 : (*capacity * 2);
        KymoKnotNativeSegment *grown = (KymoKnotNativeSegment *)realloc(
            out->segments, sizeof(KymoKnotNativeSegment) * (size_t)next);
        if (grown == NULL)
            return 0;
        out->segments = grown;
        *capacity = next;
    }

    KymoKnotNativeSegment *segment = &out->segments[out->count];
    memset(segment, 0, sizeof(*segment));
    strncpy(segment->knotType, knotType, sizeof(segment->knotType) - 1);
    segment->start = start;
    segment->end = end;
    segment->length = length;
    ++out->count;
    return 1;
}

int kymoknot_native_analyze_linear(
    const double *xyz,
    int n,
    int mode,
    unsigned int seed,
    KymoKnotNativeResult *out)
{
    if (out == NULL)
        return 0;

    out->available = 0;
    out->success = 0;
    out->count = 0;
    out->segments = NULL;

    if (xyz == NULL || n < 2)
        return 0;

    /* KymoKnot rejects chains whose endpoints coincide (they are rings). */
    {
        const double dx = xyz[0] - xyz[3 * (n - 1) + 0];
        const double dy = xyz[1] - xyz[3 * (n - 1) + 1];
        const double dz = xyz[2] - xyz[3 * (n - 1) + 2];
        if (sqrt(dx * dx + dy * dy + dz * dz) < 1.0e-6)
            return 0;
    }

    KNTID_wspace *ws = ensure_workspace(n);
    if (ws == NULL)
        return 0;

    KNTID_unset_target_knot(ws);

    {
        unsigned int mixed = seed * 2654435761u + (++g_seed_counter) * 2246822519u;
        if (mixed == 0)
            mixed = 1;
        initran_(&mixed);
    }

    KNTarc *arc_ptr = build_linear_arc(xyz, n);
    if (arc_ptr == NULL)
        return 0;
    if (arc_ptr->coord == NULL || arc_ptr->index == NULL)
    {
        KNTfree_arc(arc_ptr);
        return 0;
    }

    const int max_stride = n / 50;
    const int min_stride = 2;

    qhull_init(3 * (n + MAX_CLS_POINTS));
    KNTarc closed_arc = KNTCqhull_hybrid_close_subchain(arc_ptr, 0, n - 1);
    KNTLcentre_rotate_random(&closed_arc);

    KNTarc *knt_rect = KNTLrectify_coord_local(
        &closed_arc, 0, n - 1, min_stride, max_stride);
    if (knt_rect == NULL)
    {
        KNTfree_arc(&closed_arc);
        KNTfree_arc(arc_ptr);
        qhull_terminate();
        return 0;
    }

    KNTid knot_type = KNTID_identify_knot(knt_rect, ws);
    closed_arc.knot_type = knt_rect->knot_type;

    const int start = get_idx_rect_chain_local(knt_rect, 0);
    const int end = get_idx_rect_chain_local(knt_rect, n - 1);

    out->available = 1;
    int capacity = 0;

    if (knot_type.k_id == K_Un)
    {
        if (!append_segment(out, &capacity, "UN", -1, -1, -1))
        {
            kymoknot_native_free_result(out);
            out->available = 1;
            out->success = 0;
            KNTfree_arc(knt_rect);
            KNTfree_arc(&closed_arc);
            KNTfree_arc(arc_ptr);
            qhull_terminate();
            return 1;
        }
        out->success = 1;
    }
    else
    {
        if (mode == 1)
        {
            KNTL_BUloc_rect1(
                &closed_arc, knt_rect, start, end,
                KNTCqhull_hybrid_close_subchain, ws);
        }
        else
        {
            KNTL_TDloc_rect1(
                &closed_arc, knt_rect, start, end,
                KNTCqhull_hybrid_close_subchain, ws);
        }

        const char resultArcType =
            (mode == 1) ? (char)ARC_ID_SK : (char)ARC_ID_SCK;
        int found = 0;

        /* The historical text parser kept only the first result per line. */
        for (KNTarc *node = &closed_arc; node != NULL; node = node->next)
        {
            if (node->arc_type != resultArcType)
                continue;

            const int seg_start = node->start;
            const int seg_end = node->end;
            const int seg_length = (seg_end > seg_start)
                ? (seg_end - seg_start)
                : (seg_end - seg_start + closed_arc.len);

            char knot_ids[128];
            knot_ids[0] = '\0';
            KNTID_get_rolfsen(knot_ids, (int)sizeof(knot_ids), node->knot_type);

            if (!append_segment(
                    out, &capacity, knot_ids, seg_start, seg_end, seg_length))
            {
                kymoknot_native_free_result(out);
                out->available = 1;
                out->success = 0;
                KNTfree_arc(knt_rect);
                KNTfree_arc(&closed_arc);
                KNTfree_arc(arc_ptr);
                qhull_terminate();
                return 1;
            }
            found = 1;
            break;
        }
        out->success = found;
    }

    KNTfree_arc(&closed_arc);
    KNTfree_arc(arc_ptr);
    KNTfree_arc(knt_rect);
    qhull_terminate();
    return 1;
}

void kymoknot_native_free_result(KymoKnotNativeResult *out)
{
    if (out == NULL)
        return;
    free(out->segments);
    out->segments = NULL;
    out->count = 0;
    out->available = 0;
    out->success = 0;
}
