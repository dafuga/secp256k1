/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_BATCH_BACKEND_H
#define SECP256K1_HARBOR_BATCH_BACKEND_H

/* Architecture backends implement the same fixed-width field interface. Keep
 * selection here so recovery and group code do not depend on ISA details. */
#if defined(HARBOR_SECP256K1_BATCH_BACKEND_ARM64_NEON)
#  define HARBOR_SECP256K1_BATCH_BACKEND_ENABLED 1
#  define HARBOR_SECP256K1_BATCH_BACKEND_WIDTH 4
#  define HARBOR_SECP256K1_BATCH_BACKEND_NAME "arm64-neon"
#elif defined(HARBOR_SECP256K1_BATCH_BACKEND_PORTABLE) || \
      defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
#  define HARBOR_SECP256K1_BATCH_BACKEND_ENABLED 1
#  define HARBOR_SECP256K1_BATCH_BACKEND_WIDTH HARBOR_SECP256K1_ARM64_SQRT_BATCH_WIDTH
#  define HARBOR_SECP256K1_BATCH_BACKEND_NAME "portable-interleaved"
#endif

#endif /* SECP256K1_HARBOR_BATCH_BACKEND_H */
