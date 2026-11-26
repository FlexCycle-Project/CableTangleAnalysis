#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>

// Standard-library distributions are implementation-defined: the same
// std::mt19937 seed can produce different sequences with libstdc++, libc++,
// or different versions of either. These adapters map the standardized
// mt19937 output directly, so a seed has the same meaning everywhere.
namespace DeterministicRandom
{
    inline double unitDouble(std::mt19937& rng)
    {
        // Two 32-bit draws -> 53 mantissa bits, the usual mt19937 mapping.
        const std::uint64_t hi = static_cast<std::uint64_t>(rng()) >> 5;
        const std::uint64_t lo = static_cast<std::uint64_t>(rng()) >> 6;
        return (static_cast<double>(hi) * 67108864.0 + static_cast<double>(lo)) *
            (1.0 / 9007199254740992.0);
    }

    template <typename Real>
    inline Real unitReal(std::mt19937& rng)
    {
        return static_cast<Real>(unitDouble(rng));
    }

    inline std::size_t uniformIndex(std::mt19937& rng, std::size_t count)
    {
        if (count <= 1)
            return 0;

        const std::uint64_t range = static_cast<std::uint64_t>(count);
        const std::uint64_t limit =
            std::numeric_limits<std::uint64_t>::max() -
            (std::numeric_limits<std::uint64_t>::max() % range);

        std::uint64_t value = 0;
        do
        {
            value = (static_cast<std::uint64_t>(rng()) << 32) |
                static_cast<std::uint64_t>(rng());
        } while (value >= limit);

        return static_cast<std::size_t>(value % range);
    }
}
