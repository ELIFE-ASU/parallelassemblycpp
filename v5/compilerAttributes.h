#pragma once

#if defined(_MSC_VER)
    #define ASSEMBLYCPP_ALWAYS_INLINE __forceinline
    #define ASSEMBLYCPP_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
    #define ASSEMBLYCPP_ALWAYS_INLINE [[gnu::always_inline]]
    #define ASSEMBLYCPP_NOINLINE [[gnu::noinline]]
#else
    #define ASSEMBLYCPP_ALWAYS_INLINE
    #define ASSEMBLYCPP_NOINLINE
#endif
