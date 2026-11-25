#ifndef KYMOKNOT_NATIVE_BRIDGE_H
#define KYMOKNOT_NATIVE_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct KymoKnotNativeSegment
{
    char knotType[32];
    int start;
    int end;
    int length;
} KymoKnotNativeSegment;

typedef struct KymoKnotNativeResult
{
    int available;
    int success;
    int count;
    KymoKnotNativeSegment *segments;
} KymoKnotNativeResult;

/*
 * Analyze one linear chain in memory.  mode: 0 = top-down, 1 = bottom-up.
 * seed: base RNG seed (0 is allowed; the bridge derives a non-zero seed).
 * Returns 1 if an analysis was produced, 0 on invalid input/allocation failure.
 * The analysis itself is stochastic, as in the original command-line tool.
 */
int kymoknot_native_analyze_linear(
    const double *xyz,
    int n,
    int mode,
    unsigned int seed,
    KymoKnotNativeResult *out);

void kymoknot_native_free_result(KymoKnotNativeResult *out);

#ifdef __cplusplus
}
#endif

#endif
