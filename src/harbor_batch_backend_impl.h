/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_BATCH_BACKEND_IMPL_H
#define SECP256K1_HARBOR_BATCH_BACKEND_IMPL_H

#include "harbor_batch_backend.h"

#if defined(HARBOR_SECP256K1_BATCH_BACKEND_ARM64_NEON)
#  include "harbor_batch_backend_arm64_neon_impl.h"
#elif defined(HARBOR_SECP256K1_BATCH_BACKEND_ENABLED)
#  include "harbor_batch_backend_portable_impl.h"
#endif

#endif /* SECP256K1_HARBOR_BATCH_BACKEND_IMPL_H */
