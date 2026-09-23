#pragma once

#if defined(_MSC_VER)
    #define PARALLELASSEMBLYCPP_ALWAYS_INLINE __forceinline
    #define PARALLELASSEMBLYCPP_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    #define PARALLELASSEMBLYCPP_ALWAYS_INLINE [[gnu::always_inline]] inline
    #define PARALLELASSEMBLYCPP_NOINLINE [[gnu::noinline]]
#else
    #define PARALLELASSEMBLYCPP_ALWAYS_INLINE inline
    #define PARALLELASSEMBLYCPP_NOINLINE
#endif
