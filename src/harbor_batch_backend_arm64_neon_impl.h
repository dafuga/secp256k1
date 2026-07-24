/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_BATCH_BACKEND_ARM64_NEON_IMPL_H
#define SECP256K1_HARBOR_BATCH_BACKEND_ARM64_NEON_IMPL_H

#if !defined(__aarch64__) || !defined(SECP256K1_WIDEMUL_INT128)
#  error "The Harbor ARM64 NEON backend requires AArch64 and the 5x52 field representation"
#endif

#include <arm_neon.h>

/* Four independent field elements in base 2^26, transposed by limb. Inputs
 * cross the 5x52/10x26 boundary once per complete exponentiation chain. */
typedef struct {
    uint32x4_t n[10];
} secp256k1_harbor_fe4;

typedef struct {
    uint64x2_t lo;
    uint64x2_t hi;
} secp256k1_harbor_u64x4;

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_zero(void) {
    secp256k1_harbor_u64x4 r;
    r.lo = vdupq_n_u64(0);
    r.hi = vdupq_n_u64(0);
    return r;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_add(
        secp256k1_harbor_u64x4 a, secp256k1_harbor_u64x4 b) {
    a.lo = vaddq_u64(a.lo, b.lo);
    a.hi = vaddq_u64(a.hi, b.hi);
    return a;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_sub(
        secp256k1_harbor_u64x4 a, secp256k1_harbor_u64x4 b) {
    a.lo = vsubq_u64(a.lo, b.lo);
    a.hi = vsubq_u64(a.hi, b.hi);
    return a;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_shl(
        secp256k1_harbor_u64x4 a, const int bits) {
    int64x2_t shift = vdupq_n_s64(bits);
    a.lo = vshlq_u64(a.lo, shift);
    a.hi = vshlq_u64(a.hi, shift);
    return a;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_shr(
        secp256k1_harbor_u64x4 a, const int bits) {
    int64x2_t shift = vdupq_n_s64(-bits);
    a.lo = vshlq_u64(a.lo, shift);
    a.hi = vshlq_u64(a.hi, shift);
    return a;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_u64x4_and(
        secp256k1_harbor_u64x4 a, uint64_t mask) {
    uint64x2_t m = vdupq_n_u64(mask);
    a.lo = vandq_u64(a.lo, m);
    a.hi = vandq_u64(a.hi, m);
    return a;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_mul_u32x4(
        uint32x4_t a, uint32x4_t b) {
    secp256k1_harbor_u64x4 r;
    r.lo = vmull_u32(vget_low_u32(a), vget_low_u32(b));
    r.hi = vmull_high_u32(a, b);
    return r;
}

static SECP256K1_INLINE void secp256k1_harbor_accum_mul_u32x4(
        secp256k1_harbor_u64x4 *acc, uint32x4_t a, uint32x4_t b, int doubled) {
    secp256k1_harbor_u64x4 p = secp256k1_harbor_mul_u32x4(a, b);
    if (doubled) p = secp256k1_harbor_u64x4_shl(p, 1);
    *acc = secp256k1_harbor_u64x4_add(*acc, p);
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_mul_r0(
        secp256k1_harbor_u64x4 a) {
    secp256k1_harbor_u64x4 r;
    r.lo = vmlal_n_u32(vshlq_n_u64(vmull_n_u32(vmovn_u64(vshrq_n_u64(a.lo, 32)), 0x3d10U), 32),
                       vmovn_u64(a.lo), 0x3d10U);
    r.hi = vmlal_n_u32(vshlq_n_u64(vmull_n_u32(vmovn_u64(vshrq_n_u64(a.hi, 32)), 0x3d10U), 32),
                       vmovn_u64(a.hi), 0x3d10U);
    return r;
}

static SECP256K1_INLINE secp256k1_harbor_u64x4 secp256k1_harbor_mul_977(
        secp256k1_harbor_u64x4 a) {
    secp256k1_harbor_u64x4 r;
    r.lo = vmlal_n_u32(vshlq_n_u64(vmull_n_u32(vmovn_u64(vshrq_n_u64(a.lo, 32)), 977U), 32),
                       vmovn_u64(a.lo), 977U);
    r.hi = vmlal_n_u32(vshlq_n_u64(vmull_n_u32(vmovn_u64(vshrq_n_u64(a.hi, 32)), 977U), 32),
                       vmovn_u64(a.hi), 977U);
    return r;
}

static SECP256K1_INLINE void secp256k1_harbor_reduce_product(
        secp256k1_harbor_fe4 *r, secp256k1_harbor_u64x4 c[19]) {
    const uint64_t mask26 = 0x3ffffffULL;
    const uint64_t mask22 = 0x3fffffULL;
    int i;

    for (i = 0; i < 18; ++i) {
        secp256k1_harbor_u64x4 carry = secp256k1_harbor_u64x4_shr(c[i], 26);
        c[i] = secp256k1_harbor_u64x4_and(c[i], mask26);
        c[i + 1] = secp256k1_harbor_u64x4_add(c[i + 1], carry);
    }

    for (i = 18; i >= 10; --i) {
        c[i - 10] = secp256k1_harbor_u64x4_add(c[i - 10], secp256k1_harbor_mul_r0(c[i]));
        c[i - 9] = secp256k1_harbor_u64x4_add(c[i - 9], secp256k1_harbor_u64x4_shl(c[i], 10));
    }

    for (i = 0; i < 9; ++i) {
        secp256k1_harbor_u64x4 carry = secp256k1_harbor_u64x4_shr(c[i], 26);
        c[i] = secp256k1_harbor_u64x4_and(c[i], mask26);
        c[i + 1] = secp256k1_harbor_u64x4_add(c[i + 1], carry);
    }

    for (i = 0; i < 2; ++i) {
        secp256k1_harbor_u64x4 high = secp256k1_harbor_u64x4_shr(c[9], 22);
        int j;
        c[9] = secp256k1_harbor_u64x4_and(c[9], mask22);
        c[0] = secp256k1_harbor_u64x4_add(c[0], secp256k1_harbor_mul_977(high));
        c[1] = secp256k1_harbor_u64x4_add(c[1], secp256k1_harbor_u64x4_shl(high, 6));
        for (j = 0; j < 9; ++j) {
            secp256k1_harbor_u64x4 carry = secp256k1_harbor_u64x4_shr(c[j], 26);
            c[j] = secp256k1_harbor_u64x4_and(c[j], mask26);
            c[j + 1] = secp256k1_harbor_u64x4_add(c[j + 1], carry);
        }
    }

    for (i = 0; i < 10; ++i) {
        r->n[i] = vcombine_u32(vmovn_u64(c[i].lo), vmovn_u64(c[i].hi));
    }
}

static SECP256K1_INLINE void secp256k1_harbor_arm64_fe4_mul(
        secp256k1_harbor_fe4 *r, const secp256k1_harbor_fe4 *a,
        const secp256k1_harbor_fe4 *b) {
    secp256k1_harbor_u64x4 c[19];
    int i, j;
    for (i = 0; i < 19; ++i) c[i] = secp256k1_harbor_u64x4_zero();
    for (i = 0; i < 10; ++i) {
        for (j = 0; j < 10; ++j) {
            secp256k1_harbor_accum_mul_u32x4(&c[i + j], a->n[i], b->n[j], 0);
        }
    }
    secp256k1_harbor_reduce_product(r, c);
}

static SECP256K1_INLINE void secp256k1_harbor_arm64_fe4_sqr(
        secp256k1_harbor_fe4 *r, const secp256k1_harbor_fe4 *a) {
    secp256k1_harbor_u64x4 c[19];
    int i, j;
    for (i = 0; i < 19; ++i) c[i] = secp256k1_harbor_u64x4_zero();
    for (i = 0; i < 10; ++i) {
        secp256k1_harbor_accum_mul_u32x4(&c[2 * i], a->n[i], a->n[i], 0);
        for (j = i + 1; j < 10; ++j) {
            secp256k1_harbor_accum_mul_u32x4(&c[i + j], a->n[i], a->n[j], 1);
        }
    }
    secp256k1_harbor_reduce_product(r, c);
}

static void secp256k1_harbor_arm64_fe4_from_fe(
        secp256k1_harbor_fe4 *r,
        const secp256k1_fe a[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH]) {
    secp256k1_fe normalized[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    uint32_t lanes[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    size_t limb, lane;
    for (lane = 0; lane < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++lane) {
        normalized[lane] = a[lane];
        secp256k1_fe_normalize(&normalized[lane]);
    }
    for (limb = 0; limb < 10; ++limb) {
        for (lane = 0; lane < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++lane) {
            lanes[lane] = (uint32_t)((normalized[lane].n[limb / 2] >> (26 * (limb & 1))) & 0x3ffffffULL);
        }
        r->n[limb] = vld1q_u32(lanes);
    }
}

static void secp256k1_harbor_arm64_fe4_to_fe(
        secp256k1_fe r[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH],
        const secp256k1_harbor_fe4 *a) {
    uint32_t limbs[10][HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    size_t limb, lane;
    for (limb = 0; limb < 10; ++limb) vst1q_u32(limbs[limb], a->n[limb]);
    for (lane = 0; lane < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++lane) {
        for (limb = 0; limb < 5; ++limb) {
            r[lane].n[limb] = (uint64_t)limbs[2 * limb][lane] |
                              ((uint64_t)limbs[2 * limb + 1][lane] << 26);
        }
#ifdef VERIFY
        r[lane].magnitude = 1;
        r[lane].normalized = 0;
#endif
        secp256k1_fe_normalize(&r[lane]);
    }
}

static int secp256k1_harbor_batch_fe_sqrt(
        secp256k1_fe r[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH],
        const secp256k1_fe a[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH]) {
    secp256k1_harbor_fe4 av, x2, x3, x6, x9, x11, x22, x44, x88;
    secp256k1_harbor_fe4 x176, x220, x223, t1, check;
    secp256k1_fe check_fe[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe normalized_a[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    size_t i;
    int j;
    int ret = 1;

    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        VERIFY_CHECK(&r[i] != &a[i]);
        SECP256K1_FE_VERIFY(&a[i]);
        SECP256K1_FE_VERIFY_MAGNITUDE(&a[i], 8);
        normalized_a[i] = a[i];
        secp256k1_fe_normalize(&normalized_a[i]);
    }
    secp256k1_harbor_arm64_fe4_from_fe(&av, normalized_a);
    secp256k1_harbor_arm64_fe4_sqr(&x2, &av);
    secp256k1_harbor_arm64_fe4_mul(&x2, &x2, &av);
    secp256k1_harbor_arm64_fe4_sqr(&x3, &x2);
    secp256k1_harbor_arm64_fe4_mul(&x3, &x3, &av);
    x6 = x3;
    for (j = 0; j < 3; ++j) secp256k1_harbor_arm64_fe4_sqr(&x6, &x6);
    secp256k1_harbor_arm64_fe4_mul(&x6, &x6, &x3);
    x9 = x6;
    for (j = 0; j < 3; ++j) secp256k1_harbor_arm64_fe4_sqr(&x9, &x9);
    secp256k1_harbor_arm64_fe4_mul(&x9, &x9, &x3);
    x11 = x9;
    for (j = 0; j < 2; ++j) secp256k1_harbor_arm64_fe4_sqr(&x11, &x11);
    secp256k1_harbor_arm64_fe4_mul(&x11, &x11, &x2);
    x22 = x11;
    for (j = 0; j < 11; ++j) secp256k1_harbor_arm64_fe4_sqr(&x22, &x22);
    secp256k1_harbor_arm64_fe4_mul(&x22, &x22, &x11);
    x44 = x22;
    for (j = 0; j < 22; ++j) secp256k1_harbor_arm64_fe4_sqr(&x44, &x44);
    secp256k1_harbor_arm64_fe4_mul(&x44, &x44, &x22);
    x88 = x44;
    for (j = 0; j < 44; ++j) secp256k1_harbor_arm64_fe4_sqr(&x88, &x88);
    secp256k1_harbor_arm64_fe4_mul(&x88, &x88, &x44);
    x176 = x88;
    for (j = 0; j < 88; ++j) secp256k1_harbor_arm64_fe4_sqr(&x176, &x176);
    secp256k1_harbor_arm64_fe4_mul(&x176, &x176, &x88);
    x220 = x176;
    for (j = 0; j < 44; ++j) secp256k1_harbor_arm64_fe4_sqr(&x220, &x220);
    secp256k1_harbor_arm64_fe4_mul(&x220, &x220, &x44);
    x223 = x220;
    for (j = 0; j < 3; ++j) secp256k1_harbor_arm64_fe4_sqr(&x223, &x223);
    secp256k1_harbor_arm64_fe4_mul(&x223, &x223, &x3);
    t1 = x223;
    for (j = 0; j < 23; ++j) secp256k1_harbor_arm64_fe4_sqr(&t1, &t1);
    secp256k1_harbor_arm64_fe4_mul(&t1, &t1, &x22);
    for (j = 0; j < 6; ++j) secp256k1_harbor_arm64_fe4_sqr(&t1, &t1);
    secp256k1_harbor_arm64_fe4_mul(&t1, &t1, &x2);
    secp256k1_harbor_arm64_fe4_sqr(&t1, &t1);
    secp256k1_harbor_arm64_fe4_sqr(&t1, &t1);
    secp256k1_harbor_arm64_fe4_sqr(&check, &t1);
    secp256k1_harbor_arm64_fe4_to_fe(r, &t1);
    secp256k1_harbor_arm64_fe4_to_fe(check_fe, &check);
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        ret &= secp256k1_fe_equal(&check_fe[i], &normalized_a[i]);
    }
    return ret;
}

#endif /* SECP256K1_HARBOR_BATCH_BACKEND_ARM64_NEON_IMPL_H */
