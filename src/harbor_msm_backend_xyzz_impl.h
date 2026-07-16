/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_MSM_BACKEND_XYZZ_IMPL_H
#define SECP256K1_HARBOR_MSM_BACKEND_XYZZ_IMPL_H

/* XYZZ caches Z^2 and Z^3 for each Pippenger bucket. Mixed affine adds then
 * avoid the Z square performed by the ordinary Jacobian formula. */
static void secp256k1_harbor_xyzz_set_infinity(secp256k1_harbor_xyzz *r) {
    r->infinity = 1;
}

static void secp256k1_harbor_xyzz_set_ge(
        secp256k1_harbor_xyzz *r, const secp256k1_ge *a) {
    if (a->infinity) {
        secp256k1_harbor_xyzz_set_infinity(r);
        return;
    }
    r->x = a->x;
    r->y = a->y;
    secp256k1_fe_set_int(&r->zz, 1);
    secp256k1_fe_set_int(&r->zzz, 1);
    r->infinity = 0;
}

static void secp256k1_harbor_xyzz_set_ge_signed(
        secp256k1_harbor_xyzz *r, const secp256k1_ge *a, int negate_a) {
    secp256k1_harbor_xyzz_set_ge(r, a);
    if (negate_a && !r->infinity) {
        secp256k1_fe_normalize_weak(&r->y);
        secp256k1_fe_negate(&r->y, &r->y, 1);
    }
}

static void secp256k1_harbor_xyzz_to_gej(
        secp256k1_gej *r, const secp256k1_harbor_xyzz *a) {
    secp256k1_fe zinv, zinv2, zinv3;
    if (a->infinity) {
        secp256k1_gej_set_infinity(r);
        return;
    }

    /* z^-1 = Z^2 / Z^3. Normalize to affine once at the end of the MSM. */
    zinv = a->zzz;
    secp256k1_fe_inv_var(&zinv, &zinv);
    secp256k1_fe_mul(&zinv, &zinv, &a->zz);
    secp256k1_fe_sqr(&zinv2, &zinv);
    secp256k1_fe_mul(&zinv3, &zinv2, &zinv);
    secp256k1_fe_mul(&r->x, &a->x, &zinv2);
    secp256k1_fe_mul(&r->y, &a->y, &zinv3);
    secp256k1_fe_set_int(&r->z, 1);
    r->infinity = 0;
}

static void secp256k1_harbor_xyzz_double(
        secp256k1_harbor_xyzz *r, const secp256k1_harbor_xyzz *a) {
    secp256k1_fe l, s, s2, t, x, y, zz, zzz;
    if (a->infinity) {
        secp256k1_harbor_xyzz_set_infinity(r);
        return;
    }

    /* The same a=0 formula as secp256k1_gej_double, with cached powers:
     * Z3 = Y1*Z1, ZZ3 = Y1^2*ZZ1, ZZZ3 = Y1^3*ZZZ1. */
    secp256k1_fe_sqr(&s, &a->y);
    secp256k1_fe_mul(&zz, &a->zz, &s);
    secp256k1_fe_mul(&zzz, &s, &a->y);
    secp256k1_fe_mul(&zzz, &zzz, &a->zzz);
    secp256k1_fe_sqr(&l, &a->x);
    secp256k1_fe_mul_int(&l, 3);
    secp256k1_fe_half(&l);
    secp256k1_fe_negate(&t, &s, 1);
    secp256k1_fe_mul(&t, &t, &a->x);
    secp256k1_fe_sqr(&x, &l);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_sqr(&s2, &s);
    secp256k1_fe_add(&t, &x);
    secp256k1_fe_mul(&y, &t, &l);
    secp256k1_fe_add(&y, &s2);
    secp256k1_fe_negate(&y, &y, 2);

    r->x = x;
    r->y = y;
    r->zz = zz;
    r->zzz = zzz;
    r->infinity = 0;
}

static void secp256k1_harbor_xyzz_add_ge_signed_noninfinity(
        secp256k1_harbor_xyzz *r, const secp256k1_harbor_xyzz *a,
        const secp256k1_ge *b, int negate_b) {
    secp256k1_fe u1, u2, s1, s2, h, i, h2, h3, t;
    secp256k1_fe x, y, zz, zzz;
    u1 = a->x;
    secp256k1_fe_mul(&u2, &b->x, &a->zz);
    s1 = a->y;
    secp256k1_fe_mul(&s2, &b->y, &a->zzz);
    if (negate_b) {
        secp256k1_fe_negate(&s2, &s2, 1);
    }
    secp256k1_fe_negate(&h, &u1, SECP256K1_GEJ_X_MAGNITUDE_MAX);
    secp256k1_fe_add(&h, &u2);
    secp256k1_fe_negate(&i, &s2, 1);
    secp256k1_fe_add(&i, &s1);
    if (secp256k1_fe_normalizes_to_zero_var(&h)) {
        if (secp256k1_fe_normalizes_to_zero_var(&i)) {
            secp256k1_harbor_xyzz_double(r, a);
        } else {
            secp256k1_harbor_xyzz_set_infinity(r);
        }
        return;
    }

    secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_mul(&zz, &a->zz, &h2);
    secp256k1_fe_mul(&h3, &h2, &h);
    secp256k1_fe_mul(&zzz, &a->zzz, &h3);
    secp256k1_fe_negate(&h2, &h2, 1);
    secp256k1_fe_negate(&h3, &h3, 1);
    secp256k1_fe_mul(&t, &u1, &h2);
    secp256k1_fe_sqr(&x, &i);
    secp256k1_fe_add(&x, &h3);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_add(&t, &x);
    secp256k1_fe_mul(&y, &t, &i);
    secp256k1_fe_mul(&h3, &h3, &s1);
    secp256k1_fe_add(&y, &h3);

    r->x = x;
    r->y = y;
    r->zz = zz;
    r->zzz = zzz;
    r->infinity = 0;
}

static void secp256k1_harbor_xyzz_add_ge_signed(
        secp256k1_harbor_xyzz *r, const secp256k1_harbor_xyzz *a,
        const secp256k1_ge *b, int negate_b) {
    if (a->infinity) {
        secp256k1_harbor_xyzz_set_ge_signed(r, b, negate_b);
        return;
    }
    if (b->infinity) {
        *r = *a;
        return;
    }
    secp256k1_harbor_xyzz_add_ge_signed_noninfinity(r, a, b, negate_b);
}

static void secp256k1_harbor_xyzz_add(
        secp256k1_harbor_xyzz *r, const secp256k1_harbor_xyzz *a,
        const secp256k1_harbor_xyzz *b) {
    secp256k1_fe u1, u2, s1, s2, h, i, h2, h3, t, zz_product, zzz_product;
    secp256k1_fe x, y, zz, zzz;
    if (a->infinity) {
        *r = *b;
        return;
    }
    if (b->infinity) {
        *r = *a;
        return;
    }

    secp256k1_fe_mul(&u1, &a->x, &b->zz);
    secp256k1_fe_mul(&u2, &b->x, &a->zz);
    secp256k1_fe_mul(&s1, &a->y, &b->zzz);
    secp256k1_fe_mul(&s2, &b->y, &a->zzz);
    secp256k1_fe_negate(&h, &u1, 1);
    secp256k1_fe_add(&h, &u2);
    secp256k1_fe_negate(&i, &s2, 1);
    secp256k1_fe_add(&i, &s1);
    if (secp256k1_fe_normalizes_to_zero_var(&h)) {
        if (secp256k1_fe_normalizes_to_zero_var(&i)) {
            secp256k1_harbor_xyzz_double(r, a);
        } else {
            secp256k1_harbor_xyzz_set_infinity(r);
        }
        return;
    }

    secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_mul(&zz_product, &a->zz, &b->zz);
    secp256k1_fe_mul(&zz, &zz_product, &h2);
    secp256k1_fe_mul(&h3, &h2, &h);
    secp256k1_fe_mul(&zzz_product, &a->zzz, &b->zzz);
    secp256k1_fe_mul(&zzz, &zzz_product, &h3);
    secp256k1_fe_negate(&h2, &h2, 1);
    secp256k1_fe_negate(&h3, &h3, 1);
    secp256k1_fe_mul(&t, &u1, &h2);
    secp256k1_fe_sqr(&x, &i);
    secp256k1_fe_add(&x, &h3);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_add(&x, &t);
    secp256k1_fe_add(&t, &x);
    secp256k1_fe_mul(&y, &t, &i);
    secp256k1_fe_mul(&h3, &h3, &s1);
    secp256k1_fe_add(&y, &h3);

    r->x = x;
    r->y = y;
    r->zz = zz;
    r->zzz = zzz;
    r->infinity = 0;
}

static int secp256k1_harbor_ecmult_pippenger_wnaf(
        secp256k1_gej *unused_buckets, int bucket_window,
        struct secp256k1_pippenger_state *state, secp256k1_gej *r,
        const secp256k1_scalar *sc, const secp256k1_ge *pt, size_t num) {
    const size_t n_wnaf = WNAF_SIZE(bucket_window + 1);
    const size_t bucket_count = ECMULT_TABLE_SIZE(bucket_window + 2);
    secp256k1_harbor_xyzz result;
    size_t no = 0;
    size_t np;
    int window;
    (void)unused_buckets;

    for (np = 0; np < num; ++np) {
        if (secp256k1_scalar_is_zero(&sc[np]) || secp256k1_ge_is_infinity(&pt[np])) continue;
        state->ps[no].input_pos = np;
        state->ps[no].skew_na = secp256k1_wnaf_fixed(
            &state->wnaf_na[no * n_wnaf], &sc[np], bucket_window + 1);
        ++no;
    }
    secp256k1_harbor_xyzz_set_infinity(&result);
    if (no == 0) {
        secp256k1_gej_set_infinity(r);
        return 1;
    }

    for (window = (int)n_wnaf - 1; window >= 0; --window) {
        secp256k1_harbor_xyzz running_sum;
        size_t bucket;
        int j;

        for (bucket = 0; bucket < bucket_count; ++bucket) {
            secp256k1_harbor_xyzz_set_infinity(&state->xyzz_buckets[bucket]);
        }
        for (np = 0; np < no; ++np) {
            const struct secp256k1_pippenger_point_state point_state = state->ps[np];
            const int digit = state->wnaf_na[np * n_wnaf + window];
#if defined(HARBOR_SECP256K1_XYZZ_PREFETCH_DISTANCE) && HARBOR_SECP256K1_XYZZ_PREFETCH_DISTANCE > 0 && defined(__GNUC__)
            if (np + HARBOR_SECP256K1_XYZZ_PREFETCH_DISTANCE < no) {
                const size_t prefetch_np = np + HARBOR_SECP256K1_XYZZ_PREFETCH_DISTANCE;
                const int prefetch_digit = state->wnaf_na[prefetch_np * n_wnaf + window];
                if (prefetch_digit != 0) {
                    const size_t prefetch_bucket = prefetch_digit > 0
                        ? (size_t)(prefetch_digit - 1) / 2
                        : (size_t)(-(prefetch_digit + 1)) / 2;
                    __builtin_prefetch(&state->xyzz_buckets[prefetch_bucket], 1, 3);
                    __builtin_prefetch((const unsigned char *)&state->xyzz_buckets[prefetch_bucket] + 128,
                                       1, 3);
                }
            }
#endif
            if (window == 0 && point_state.skew_na) {
                secp256k1_harbor_xyzz_add_ge_signed(
                    &state->xyzz_buckets[0], &state->xyzz_buckets[0],
                    &pt[point_state.input_pos], 1);
            }
            if (digit == 0) continue;
            if (digit > 0) {
                bucket = (size_t)(digit - 1) / 2;
            } else {
                bucket = (size_t)(-(digit + 1)) / 2;
            }
            secp256k1_harbor_xyzz_add_ge_signed(
                &state->xyzz_buckets[bucket], &state->xyzz_buckets[bucket],
                &pt[point_state.input_pos], digit < 0);
        }

        for (j = 0; j < bucket_window; ++j) {
            secp256k1_harbor_xyzz_double(&result, &result);
        }
        secp256k1_harbor_xyzz_set_infinity(&running_sum);
        for (bucket = bucket_count - 1; bucket > 0; --bucket) {
            secp256k1_harbor_xyzz_add(
                &running_sum, &running_sum, &state->xyzz_buckets[bucket]);
            secp256k1_harbor_xyzz_add(&result, &result, &running_sum);
        }
        secp256k1_harbor_xyzz_add(
            &running_sum, &running_sum, &state->xyzz_buckets[0]);
        secp256k1_harbor_xyzz_double(&result, &result);
        secp256k1_harbor_xyzz_add(&result, &result, &running_sum);
    }

    secp256k1_harbor_xyzz_to_gej(r, &result);
    return 1;
}

#endif /* SECP256K1_HARBOR_MSM_BACKEND_XYZZ_IMPL_H */
