/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_MSM_BACKEND_BATCH_AFFINE_IMPL_H
#define SECP256K1_HARBOR_MSM_BACKEND_BATCH_AFFINE_IMPL_H

static SECP256K1_INLINE void secp256k1_harbor_fe_sub(
        secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b) {
    secp256k1_fe neg = *b;
    secp256k1_fe_normalize_weak(&neg);
    secp256k1_fe_negate(&neg, &neg, 1);
    *r = *a;
    secp256k1_fe_add(r, &neg);
    secp256k1_fe_normalize_weak(r);
}

static SECP256K1_INLINE void secp256k1_harbor_affine_add_with_inverse(
        secp256k1_ge *r, const secp256k1_ge *a, const secp256k1_ge *b,
        const secp256k1_fe *inverse_dx) {
    secp256k1_fe lambda, sum_x, delta;

    secp256k1_harbor_fe_sub(&lambda, &b->y, &a->y);
    secp256k1_fe_mul(&lambda, &lambda, inverse_dx);

    sum_x = a->x;
    secp256k1_fe_add(&sum_x, &b->x);
    secp256k1_fe_normalize_weak(&sum_x);
    secp256k1_fe_negate(&sum_x, &sum_x, 1);
    secp256k1_fe_sqr(&r->x, &lambda);
    secp256k1_fe_add(&r->x, &sum_x);
    secp256k1_fe_normalize_weak(&r->x);

    secp256k1_harbor_fe_sub(&delta, &a->x, &r->x);
    secp256k1_fe_mul(&r->y, &lambda, &delta);
    secp256k1_harbor_fe_sub(&r->y, &r->y, &a->y);
    r->infinity = 0;
}

/* Returns zero when an exceptional affine pair requires the upstream complete
 * group formula. The caller then reruns the unchanged upstream MSM. */
static int secp256k1_harbor_ecmult_pippenger_wnaf(
        secp256k1_gej *buckets, int bucket_window,
        struct secp256k1_pippenger_state *state, secp256k1_gej *r,
        const secp256k1_scalar *sc, const secp256k1_ge *pt, size_t num) {
    const size_t n_wnaf = WNAF_SIZE(bucket_window + 1);
    const size_t bucket_count = ECMULT_TABLE_SIZE(bucket_window + 2);
    size_t no = 0;
    size_t np;
    int window;

    for (np = 0; np < num; ++np) {
        if (secp256k1_scalar_is_zero(&sc[np]) || secp256k1_ge_is_infinity(&pt[np])) continue;
        state->ps[no].input_pos = np;
        state->ps[no].skew_na = secp256k1_wnaf_fixed(
            &state->wnaf_na[no * n_wnaf], &sc[np], bucket_window + 1);
        ++no;
    }
    secp256k1_gej_set_infinity(r);
    if (no == 0) return 1;

    for (window = (int)n_wnaf - 1; window >= 0; --window) {
        size_t total = 0;
        size_t bucket;
        int j;

        memset(state->bucket_counts, 0, bucket_count * sizeof(*state->bucket_counts));
        for (np = 0; np < no; ++np) {
            const int digit = state->wnaf_na[np * n_wnaf + window];
            if (window == 0 && state->ps[np].skew_na) ++state->bucket_counts[0];
            if (digit > 0) {
                ++state->bucket_counts[(size_t)(digit - 1) / 2];
            } else if (digit < 0) {
                ++state->bucket_counts[(size_t)(-(digit + 1)) / 2];
            }
        }
        for (bucket = 0; bucket < bucket_count; ++bucket) {
            state->bucket_offsets[bucket] = total;
            state->bucket_write[bucket] = total;
            total += state->bucket_counts[bucket];
        }
        VERIFY_CHECK(total <= state->affine_capacity);

        for (np = 0; np < no; ++np) {
            const struct secp256k1_pippenger_point_state point_state = state->ps[np];
            const int digit = state->wnaf_na[np * n_wnaf + window];
            secp256k1_ge point;
            size_t target;
            if (window == 0 && point_state.skew_na) {
                secp256k1_ge_neg(&point, &pt[point_state.input_pos]);
                target = state->bucket_write[0]++;
                state->affine_points[target] = point;
            }
            if (digit == 0) continue;
            if (digit > 0) {
                bucket = (size_t)(digit - 1) / 2;
                point = pt[point_state.input_pos];
            } else {
                bucket = (size_t)(-(digit + 1)) / 2;
                secp256k1_ge_neg(&point, &pt[point_state.input_pos]);
            }
            target = state->bucket_write[bucket]++;
            state->affine_points[target] = point;
        }

        for (;;) {
            secp256k1_fe accumulator;
            size_t pair_count = 0;
            size_t pair_index = 0;

            for (bucket = 0; bucket < bucket_count; ++bucket) {
                pair_count += state->bucket_counts[bucket] / 2;
            }
            if (pair_count == 0) break;
            VERIFY_CHECK(pair_count <= state->inverse_capacity);

            for (bucket = 0; bucket < bucket_count; ++bucket) {
                const size_t start = state->bucket_offsets[bucket];
                const size_t pairs = state->bucket_counts[bucket] / 2;
                size_t pair;
                for (pair = 0; pair < pairs; ++pair) {
                    const secp256k1_ge *a = &state->affine_points[start + 2 * pair];
                    const secp256k1_ge *b = &state->affine_points[start + 2 * pair + 1];
                    secp256k1_harbor_fe_sub(&state->denominators[pair_index], &b->x, &a->x);
                    secp256k1_fe_normalize_var(&state->denominators[pair_index]);
                    if (secp256k1_fe_is_zero(&state->denominators[pair_index])) return 0;
                    ++pair_index;
                }
            }

            secp256k1_fe_set_int(&accumulator, 1);
            for (pair_index = 0; pair_index < pair_count; ++pair_index) {
                state->inverse_prefix[pair_index] = accumulator;
                secp256k1_fe_mul(&accumulator, &accumulator, &state->denominators[pair_index]);
            }
            secp256k1_fe_inv_var(&accumulator, &accumulator);
            pair_index = pair_count;
            while (pair_index > 0) {
                secp256k1_fe inverse;
                --pair_index;
                secp256k1_fe_mul(&inverse, &accumulator, &state->inverse_prefix[pair_index]);
                secp256k1_fe_mul(&accumulator, &accumulator, &state->denominators[pair_index]);
                state->inverse_prefix[pair_index] = inverse;
            }

            pair_index = 0;
            for (bucket = 0; bucket < bucket_count; ++bucket) {
                const size_t start = state->bucket_offsets[bucket];
                const size_t count = state->bucket_counts[bucket];
                const size_t pairs = count / 2;
                size_t pair;
                for (pair = 0; pair < pairs; ++pair) {
                    const secp256k1_ge a = state->affine_points[start + 2 * pair];
                    const secp256k1_ge b = state->affine_points[start + 2 * pair + 1];
                    secp256k1_harbor_affine_add_with_inverse(
                        &state->affine_points[start + pair], &a, &b,
                        &state->inverse_prefix[pair_index++]);
                }
                if (count & 1U) {
                    state->affine_points[start + pairs] = state->affine_points[start + count - 1];
                }
                state->bucket_counts[bucket] = pairs + (count & 1U);
            }
        }

        for (bucket = 0; bucket < bucket_count; ++bucket) {
            if (state->bucket_counts[bucket] == 0) {
                secp256k1_gej_set_infinity(&buckets[bucket]);
            } else {
                secp256k1_gej_set_ge(&buckets[bucket],
                                     &state->affine_points[state->bucket_offsets[bucket]]);
            }
        }

        for (j = 0; j < bucket_window; ++j) secp256k1_gej_double_var(r, r, NULL);
        {
            secp256k1_gej running_sum;
            secp256k1_gej_set_infinity(&running_sum);
            for (bucket = bucket_count - 1; bucket > 0; --bucket) {
                secp256k1_gej_add_var(&running_sum, &running_sum, &buckets[bucket], NULL);
                secp256k1_gej_add_var(r, r, &running_sum, NULL);
            }
            secp256k1_gej_add_var(&running_sum, &running_sum, &buckets[0], NULL);
            secp256k1_gej_double_var(r, r, NULL);
            secp256k1_gej_add_var(r, r, &running_sum, NULL);
        }
    }
    return 1;
}

#endif /* SECP256K1_HARBOR_MSM_BACKEND_BATCH_AFFINE_IMPL_H */
