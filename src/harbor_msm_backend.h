/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_MSM_BACKEND_H
#define SECP256K1_HARBOR_MSM_BACKEND_H

/* MSM policy is independent from the field ISA backend. Architecture-specific
 * kernels can reuse the same batch-affine work decomposition. */
#if defined(HARBOR_SECP256K1_MSM_BACKEND_BATCH_AFFINE)
#  define HARBOR_SECP256K1_MSM_BACKEND_ENABLED 1
#  define HARBOR_SECP256K1_MSM_BACKEND_NAME "batch-affine"
#endif

#endif /* SECP256K1_HARBOR_MSM_BACKEND_H */
