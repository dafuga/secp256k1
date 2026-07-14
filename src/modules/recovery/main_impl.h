/***********************************************************************
 * Copyright (c) 2013-2015 Pieter Wuille                               *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_MODULE_RECOVERY_MAIN_H
#define SECP256K1_MODULE_RECOVERY_MAIN_H

#include "../../../include/secp256k1_recovery.h"

static void secp256k1_ecdsa_recoverable_signature_load(const secp256k1_context* ctx, secp256k1_scalar* r, secp256k1_scalar* s, int* recid, const secp256k1_ecdsa_recoverable_signature* sig) {
    (void)ctx;
    if (sizeof(secp256k1_scalar) == 32) {
        /* When the secp256k1_scalar type is exactly 32 byte, use its
         * representation inside secp256k1_ecdsa_signature, as conversion is very fast.
         * Note that secp256k1_ecdsa_signature_save must use the same representation. */
        memcpy(r, &sig->data[0], 32);
        memcpy(s, &sig->data[32], 32);
    } else {
        secp256k1_scalar_set_b32(r, &sig->data[0], NULL);
        secp256k1_scalar_set_b32(s, &sig->data[32], NULL);
    }
    *recid = sig->data[64];
}

static void secp256k1_ecdsa_recoverable_signature_save(secp256k1_ecdsa_recoverable_signature* sig, const secp256k1_scalar* r, const secp256k1_scalar* s, int recid) {
    if (sizeof(secp256k1_scalar) == 32) {
        memcpy(&sig->data[0], r, 32);
        memcpy(&sig->data[32], s, 32);
    } else {
        secp256k1_scalar_get_b32(&sig->data[0], r);
        secp256k1_scalar_get_b32(&sig->data[32], s);
    }
    sig->data[64] = recid;
}

int secp256k1_ecdsa_recoverable_signature_parse_compact(const secp256k1_context* ctx, secp256k1_ecdsa_recoverable_signature* sig, const unsigned char *input64, int recid) {
    secp256k1_scalar r, s;
    int ret = 1;
    int overflow = 0;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(sig != NULL);
    ARG_CHECK(input64 != NULL);
    ARG_CHECK(recid >= 0 && recid <= 3);

    secp256k1_scalar_set_b32(&r, &input64[0], &overflow);
    ret &= !overflow;
    secp256k1_scalar_set_b32(&s, &input64[32], &overflow);
    ret &= !overflow;
    if (ret) {
        secp256k1_ecdsa_recoverable_signature_save(sig, &r, &s, recid);
    } else {
        memset(sig, 0, sizeof(*sig));
    }
    return ret;
}

int secp256k1_ecdsa_recoverable_signature_serialize_compact(const secp256k1_context* ctx, unsigned char *output64, int *recid, const secp256k1_ecdsa_recoverable_signature* sig) {
    secp256k1_scalar r, s;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(output64 != NULL);
    ARG_CHECK(sig != NULL);
    ARG_CHECK(recid != NULL);

    secp256k1_ecdsa_recoverable_signature_load(ctx, &r, &s, recid, sig);
    secp256k1_scalar_get_b32(&output64[0], &r);
    secp256k1_scalar_get_b32(&output64[32], &s);
    return 1;
}

int secp256k1_ecdsa_recoverable_signature_convert(const secp256k1_context* ctx, secp256k1_ecdsa_signature* sig, const secp256k1_ecdsa_recoverable_signature* sigin) {
    secp256k1_scalar r, s;
    int recid;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(sig != NULL);
    ARG_CHECK(sigin != NULL);

    secp256k1_ecdsa_recoverable_signature_load(ctx, &r, &s, &recid, sigin);
    secp256k1_ecdsa_signature_save(sig, &r, &s);
    return 1;
}

static int secp256k1_ecdsa_sig_recover(const secp256k1_scalar *sigr, const secp256k1_scalar* sigs, secp256k1_ge *pubkey, const secp256k1_scalar *message, int recid) {
    unsigned char brx[32];
    secp256k1_fe fx;
    secp256k1_ge x;
    secp256k1_gej xj;
    secp256k1_scalar rn, u1, u2;
    secp256k1_gej qj;
    int r;

    if (secp256k1_scalar_is_zero(sigr) || secp256k1_scalar_is_zero(sigs)) {
        return 0;
    }

    secp256k1_scalar_get_b32(brx, sigr);
    r = secp256k1_fe_set_b32_limit(&fx, brx);
    (void)r;
    VERIFY_CHECK(r); /* brx comes from a scalar, so is less than the order; certainly less than p */
    if (recid & 2) {
        if (secp256k1_fe_cmp_var(&fx, &secp256k1_ecdsa_const_p_minus_order) >= 0) {
            return 0;
        }
        secp256k1_fe_add(&fx, &secp256k1_ecdsa_const_order_as_fe);
    }
    if (!secp256k1_ge_set_xo_var(&x, &fx, recid & 1)) {
        return 0;
    }
    secp256k1_gej_set_ge(&xj, &x);
    secp256k1_scalar_inverse_var(&rn, sigr);
    secp256k1_scalar_mul(&u1, &rn, message);
    secp256k1_scalar_negate(&u1, &u1);
    secp256k1_scalar_mul(&u2, &rn, sigs);
    secp256k1_ecmult(&qj, &xj, &u2, &u1);
    secp256k1_ge_set_gej_var(pubkey, &qj);
    return !secp256k1_gej_is_infinity(&qj);
}

int secp256k1_ecdsa_sign_recoverable(const secp256k1_context* ctx, secp256k1_ecdsa_recoverable_signature *signature, const unsigned char *msghash32, const unsigned char *seckey, secp256k1_nonce_function noncefp, const void* noncedata) {
    secp256k1_scalar r, s;
    int ret, recid;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(secp256k1_ecmult_gen_context_is_built(&ctx->ecmult_gen_ctx));
    ARG_CHECK(msghash32 != NULL);
    ARG_CHECK(signature != NULL);
    ARG_CHECK(seckey != NULL);

    ret = secp256k1_ecdsa_sign_inner(ctx, &r, &s, &recid, msghash32, seckey, noncefp, noncedata);
    secp256k1_ecdsa_recoverable_signature_save(signature, &r, &s, recid);
    return ret;
}

int secp256k1_ecdsa_recover(const secp256k1_context* ctx, secp256k1_pubkey *pubkey, const secp256k1_ecdsa_recoverable_signature *signature, const unsigned char *msghash32) {
    secp256k1_ge q;
    secp256k1_scalar r, s;
    secp256k1_scalar m;
    int recid;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(msghash32 != NULL);
    ARG_CHECK(signature != NULL);
    ARG_CHECK(pubkey != NULL);

    secp256k1_ecdsa_recoverable_signature_load(ctx, &r, &s, &recid, signature);
    VERIFY_CHECK(recid >= 0 && recid < 4);  /* should have been caught in parse_compact */
    secp256k1_scalar_set_b32(&m, msghash32, NULL);
    if (secp256k1_ecdsa_sig_recover(&r, &s, &q, &m, recid)) {
        secp256k1_pubkey_save(pubkey, &q);
        return 1;
    } else {
        memset(pubkey, 0, sizeof(*pubkey));
        return 0;
    }
}

typedef struct {
    secp256k1_scalar *scalars;
    secp256k1_ge *points;
} secp256k1_ecdsa_batch_data;

struct secp256k1_ecdsa_recoverable_batch_workspace_struct {
    secp256k1_scalar *scalars;
    secp256k1_ge *points;
    size_t *pubkey_slots;
    size_t *pubkey_groups;
    size_t *pubkey_group_first;
    size_t pubkey_slot_capacity;
#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
    unsigned char *recid_odds;
#endif
    secp256k1_scratch *scratch;
    size_t capacity;
#ifdef HARBOR_SECP256K1_BATCH_COEFFICIENTS_128
    unsigned char coefficient_pair[32];
#endif
};

static int secp256k1_ecdsa_batch_callback(secp256k1_scalar *scalar, secp256k1_ge *point,
                                           size_t idx, void *data) {
    secp256k1_ecdsa_batch_data *batch = (secp256k1_ecdsa_batch_data *)data;
    *scalar = batch->scalars[idx];
    *point = batch->points[idx];
    return 1;
}

secp256k1_ecdsa_recoverable_batch_workspace *secp256k1_ecdsa_recoverable_batch_workspace_create(
        const secp256k1_context *ctx, size_t capacity) {
    secp256k1_ecdsa_recoverable_batch_workspace *workspace = NULL;
    size_t scratch_size;

    VERIFY_CHECK(ctx != NULL);
    if (capacity == 0 || capacity > SIZE_MAX / (2 * sizeof(*workspace->scalars)) ||
        capacity > SIZE_MAX / (2 * sizeof(*workspace->points)) ||
        capacity > (SIZE_MAX - 1024 * 1024) / 4096) {
        return NULL;
    }

    workspace = (secp256k1_ecdsa_recoverable_batch_workspace *)malloc(sizeof(*workspace));
    if (workspace == NULL) return NULL;
    memset(workspace, 0, sizeof(*workspace));
    workspace->capacity = capacity;
    workspace->scalars = (secp256k1_scalar *)malloc(2 * capacity * sizeof(*workspace->scalars));
    workspace->points = (secp256k1_ge *)malloc(2 * capacity * sizeof(*workspace->points));
    workspace->pubkey_slot_capacity = 1;
    while (workspace->pubkey_slot_capacity < 2 * capacity) {
        if (workspace->pubkey_slot_capacity > SIZE_MAX / 2) {
            secp256k1_ecdsa_recoverable_batch_workspace_destroy(ctx, workspace);
            return NULL;
        }
        workspace->pubkey_slot_capacity *= 2;
    }
    workspace->pubkey_slots = (size_t *)malloc(workspace->pubkey_slot_capacity * sizeof(*workspace->pubkey_slots));
    workspace->pubkey_groups = (size_t *)malloc(capacity * sizeof(*workspace->pubkey_groups));
    workspace->pubkey_group_first = (size_t *)malloc(capacity * sizeof(*workspace->pubkey_group_first));
#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
    workspace->recid_odds = (unsigned char *)malloc(capacity * sizeof(*workspace->recid_odds));
#endif
    scratch_size = 1024 * 1024 + capacity * 4096;
    workspace->scratch = secp256k1_scratch_create(&ctx->error_callback, scratch_size);
    if (workspace->scalars == NULL || workspace->points == NULL ||
        workspace->pubkey_slots == NULL || workspace->pubkey_groups == NULL ||
        workspace->pubkey_group_first == NULL || workspace->scratch == NULL
#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
        || workspace->recid_odds == NULL
#endif
        ) {
        secp256k1_ecdsa_recoverable_batch_workspace_destroy(ctx, workspace);
        return NULL;
    }
    return workspace;
}

void secp256k1_ecdsa_recoverable_batch_workspace_destroy(
        const secp256k1_context *ctx,
        secp256k1_ecdsa_recoverable_batch_workspace *workspace) {
    VERIFY_CHECK(ctx != NULL);
    if (workspace == NULL) return;
    if (workspace->scratch != NULL) {
        secp256k1_scratch_destroy(&ctx->error_callback, workspace->scratch);
    }
    free(workspace->scalars);
    free(workspace->points);
    free(workspace->pubkey_slots);
    free(workspace->pubkey_groups);
    free(workspace->pubkey_group_first);
#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
    free(workspace->recid_odds);
#endif
    free(workspace);
}

static uint64_t secp256k1_ecdsa_batch_pubkey_hash(const secp256k1_pubkey *pubkey) {
    const unsigned char *bytes = (const unsigned char *)pubkey;
    uint64_t hash = UINT64_C(1469598103934665603);
    size_t i;
    for (i = 0; i < sizeof(*pubkey); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    hash ^= hash >> 33;
    hash *= UINT64_C(0xff51afd7ed558ccd);
    hash ^= hash >> 33;
    return hash;
}

/* Return the number of distinct public-key terms worth aggregating. An
 * identity result means the caller should retain the original 2N equation.
 * Probe length is bounded so adversarial keys cannot turn grouping into an
 * unbounded admission cost. */
static size_t secp256k1_ecdsa_batch_group_pubkeys(
        secp256k1_ecdsa_recoverable_batch_workspace *workspace,
        const secp256k1_pubkey *pubkeys,
        size_t count) {
    const size_t mask = workspace->pubkey_slot_capacity - 1;
    size_t unique_count = 0;
    size_t i;

    if (count < 256) return count;
    memset(workspace->pubkey_slots, 0xff,
           workspace->pubkey_slot_capacity * sizeof(*workspace->pubkey_slots));

    for (i = 0; i < count; ++i) {
        size_t slot = (size_t)secp256k1_ecdsa_batch_pubkey_hash(&pubkeys[i]) & mask;
        size_t probes;
        for (probes = 0; probes < 64; ++probes) {
            const size_t first = workspace->pubkey_slots[slot];
            if (first == SIZE_MAX) {
                workspace->pubkey_slots[slot] = i;
                workspace->pubkey_groups[i] = unique_count;
                workspace->pubkey_group_first[unique_count] = i;
                ++unique_count;
                break;
            }
            if (memcmp(&pubkeys[first], &pubkeys[i], sizeof(pubkeys[i])) == 0) {
                workspace->pubkey_groups[i] = workspace->pubkey_groups[first];
                break;
            }
            slot = (slot + 1) & mask;
        }
        if (probes == 64) return count;
    }

    /* The hash pass is cheap, but compacting a nearly unique batch is not. */
    if (unique_count * 8 > count * 7) return count;
    return unique_count;
}

#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
static int secp256k1_ecdsa_batch_prepare_r_points(
        const secp256k1_context *ctx,
        secp256k1_ecdsa_recoverable_batch_workspace *workspace,
        const secp256k1_ecdsa_recoverable_signature *signatures,
        size_t count) {
    size_t i;

    for (i = 0; i < count; ++i) {
        secp256k1_scalar r, s;
        unsigned char r_bytes[32];
        int recid;

        secp256k1_ecdsa_recoverable_signature_load(ctx, &r, &s, &recid, &signatures[i]);
        if (secp256k1_scalar_is_zero(&r) || secp256k1_scalar_is_zero(&s)) return 0;
        if (secp256k1_scalar_is_high(&s)) recid ^= 1;
        secp256k1_scalar_get_b32(r_bytes, &r);
        if (!secp256k1_fe_set_b32_limit(&workspace->points[i].x, r_bytes)) return 0;
        if (recid & 2) {
            if (secp256k1_fe_cmp_var(&workspace->points[i].x,
                                     &secp256k1_ecdsa_const_p_minus_order) >= 0) return 0;
            secp256k1_fe_add(&workspace->points[i].x, &secp256k1_ecdsa_const_order_as_fe);
        }
        workspace->recid_odds[i] = (unsigned char)(recid & 1);
    }

    i = 0;
    while (count - i >= HARBOR_SECP256K1_ARM64_SQRT_BATCH_WIDTH) {
        if (!secp256k1_ge_set_xo_var_interleaved(&workspace->points[i],
                                                  &workspace->recid_odds[i])) return 0;
        i += HARBOR_SECP256K1_ARM64_SQRT_BATCH_WIDTH;
    }
    for (; i < count; ++i) {
        if (!secp256k1_ge_set_xo_var(&workspace->points[i], &workspace->points[i].x,
                                      workspace->recid_odds[i])) return 0;
    }
    return 1;
}
#endif

int secp256k1_ecdsa_recoverable_verify_batch_workspace(
                                              const secp256k1_context *ctx,
                                              secp256k1_ecdsa_recoverable_batch_workspace *workspace,
                                              const secp256k1_ecdsa_recoverable_signature *signatures,
                                              const unsigned char *messages32,
                                              const secp256k1_pubkey *pubkeys,
                                              size_t count) {
#ifdef HARBOR_SECP256K1_BATCH_COEFFICIENTS_128
    static const unsigned char domain[] = "Harbor/K1BatchVerify/128/v1";
#else
    static const unsigned char domain[] = "Harbor/K1BatchVerify/v1";
#endif
    secp256k1_ecdsa_batch_data batch;
    secp256k1_scalar g_scalar = secp256k1_scalar_zero;
    secp256k1_gej result;
    secp256k1_sha256 transcript;
#ifndef HARBOR_SECP256K1_BATCH_COEFFICIENTS_128
    secp256k1_sha256 coefficient_prefix;
#endif
    unsigned char transcript_hash[32];
    size_t grouped_pubkeys;
    int aggregate_pubkeys;
    size_t i;
    int ok = 0;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(workspace != NULL);
    ARG_CHECK(signatures != NULL);
    ARG_CHECK(messages32 != NULL);
    ARG_CHECK(pubkeys != NULL);
    if (count == 0) return 1;
    if (count > workspace->capacity) return 0;

    batch.scalars = workspace->scalars;
    batch.points = workspace->points;

    grouped_pubkeys = secp256k1_ecdsa_batch_group_pubkeys(workspace, pubkeys, count);
    aggregate_pubkeys = grouped_pubkeys < count;
    if (aggregate_pubkeys) {
        for (i = 0; i < grouped_pubkeys; ++i) {
            batch.scalars[count + i] = secp256k1_scalar_zero;
        }
    }

    secp256k1_sha256_initialize(&transcript);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, domain, sizeof(domain) - 1);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, (const unsigned char *)signatures, count * sizeof(*signatures));
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, messages32, count * 32);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, (const unsigned char *)pubkeys, count * sizeof(*pubkeys));
    secp256k1_sha256_finalize(secp256k1_get_hash_context(ctx), &transcript, transcript_hash);

#ifndef HARBOR_SECP256K1_BATCH_COEFFICIENTS_128
    secp256k1_sha256_initialize(&coefficient_prefix);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_prefix,
                           transcript_hash, sizeof(transcript_hash));
#endif

#if defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
    if (!secp256k1_ecdsa_batch_prepare_r_points(ctx, workspace, signatures, count)) goto cleanup;
    if (!aggregate_pubkeys) {
        i = count;
        while (i > 0) {
            --i;
            batch.points[2 * i] = batch.points[i];
        }
    }
#endif

    for (i = 0; i < count; ++i) {
        secp256k1_scalar r, s, z, coefficient, term, q_term;
#if !defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
        secp256k1_fe x;
#endif
        secp256k1_sha256 coefficient_hash;
        unsigned char coefficient_bytes[32];
#if !defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
        unsigned char r_bytes[32];
#endif
        unsigned char index_bytes[8];
        int recid;
        int overflow;
        size_t j;

        secp256k1_ecdsa_recoverable_signature_load(ctx, &r, &s, &recid, &signatures[i]);
        if (secp256k1_scalar_is_zero(&r) || secp256k1_scalar_is_zero(&s)) goto cleanup;
        if (secp256k1_scalar_is_high(&s)) {
            secp256k1_scalar_negate(&s, &s);
            recid ^= 1;
        }

#if !defined(HARBOR_SECP256K1_ARM64_INTERLEAVED_SQRT)
        secp256k1_scalar_get_b32(r_bytes, &r);
        if (!secp256k1_fe_set_b32_limit(&x, r_bytes)) goto cleanup;
        if (recid & 2) {
            if (secp256k1_fe_cmp_var(&x, &secp256k1_ecdsa_const_p_minus_order) >= 0) goto cleanup;
            secp256k1_fe_add(&x, &secp256k1_ecdsa_const_order_as_fe);
        }
        if (aggregate_pubkeys) {
            const size_t group = workspace->pubkey_groups[i];
            if (!secp256k1_ge_set_xo_var(&batch.points[i], &x, recid & 1)) goto cleanup;
            if (workspace->pubkey_group_first[group] == i &&
                !secp256k1_pubkey_load(ctx, &batch.points[count + group], &pubkeys[i])) goto cleanup;
        } else {
            if (!secp256k1_ge_set_xo_var(&batch.points[2 * i], &x, recid & 1)) goto cleanup;
            if (!secp256k1_pubkey_load(ctx, &batch.points[2 * i + 1], &pubkeys[i])) goto cleanup;
        }
#else
        if (aggregate_pubkeys) {
            const size_t group = workspace->pubkey_groups[i];
            if (workspace->pubkey_group_first[group] == i &&
                !secp256k1_pubkey_load(ctx, &batch.points[count + group], &pubkeys[i])) goto cleanup;
        } else if (!secp256k1_pubkey_load(ctx, &batch.points[2 * i + 1], &pubkeys[i])) {
            goto cleanup;
        }
#endif

#ifdef HARBOR_SECP256K1_BATCH_COEFFICIENTS_128
        if ((i & 1U) == 0) {
            const size_t coefficient_pair = i / 2;
            for (j = 0; j < sizeof(index_bytes); ++j)
                index_bytes[sizeof(index_bytes) - 1 - j] = (unsigned char)(coefficient_pair >> (8 * j));
            secp256k1_sha256_initialize(&coefficient_hash);
            secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_hash,
                                   transcript_hash, sizeof(transcript_hash));
            secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_hash,
                                   index_bytes, sizeof(index_bytes));
            secp256k1_sha256_finalize(secp256k1_get_hash_context(ctx), &coefficient_hash,
                                     workspace->coefficient_pair);
        }
        memset(coefficient_bytes, 0, 16);
        memcpy(coefficient_bytes + 16, workspace->coefficient_pair + ((i & 1U) ? 16 : 0), 16);
#else
        for (j = 0; j < sizeof(index_bytes); ++j)
            index_bytes[sizeof(index_bytes) - 1 - j] = (unsigned char)(i >> (8 * j));
        coefficient_hash = coefficient_prefix;
        secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_hash,
                               index_bytes, sizeof(index_bytes));
        secp256k1_sha256_finalize(secp256k1_get_hash_context(ctx), &coefficient_hash,
                                 coefficient_bytes);
#endif
        secp256k1_scalar_set_b32(&coefficient, coefficient_bytes, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&coefficient))
            secp256k1_scalar_set_int(&coefficient, 1);

        if (aggregate_pubkeys) {
            const size_t group = workspace->pubkey_groups[i];
            secp256k1_scalar_mul(&batch.scalars[i], &coefficient, &s);
            secp256k1_scalar_mul(&q_term, &coefficient, &r);
            secp256k1_scalar_negate(&q_term, &q_term);
            secp256k1_scalar_add(&batch.scalars[count + group],
                                 &batch.scalars[count + group], &q_term);
        } else {
            secp256k1_scalar_mul(&batch.scalars[2 * i], &coefficient, &s);
            secp256k1_scalar_mul(&batch.scalars[2 * i + 1], &coefficient, &r);
            secp256k1_scalar_negate(&batch.scalars[2 * i + 1], &batch.scalars[2 * i + 1]);
        }
        secp256k1_scalar_set_b32(&z, messages32 + i * 32, NULL);
        secp256k1_scalar_mul(&term, &coefficient, &z);
        secp256k1_scalar_add(&g_scalar, &g_scalar, &term);
    }
    secp256k1_scalar_negate(&g_scalar, &g_scalar);

    if (!secp256k1_ecmult_multi_var(&ctx->error_callback, workspace->scratch, &result, &g_scalar,
                                    secp256k1_ecdsa_batch_callback, &batch,
                                    aggregate_pubkeys ? count + grouped_pubkeys : 2 * count)) goto cleanup;
    ok = secp256k1_gej_is_infinity(&result);

cleanup:
    return ok;
}

int secp256k1_ecdsa_recoverable_verify_batch(const secp256k1_context *ctx,
                                              const secp256k1_ecdsa_recoverable_signature *signatures,
                                              const unsigned char *messages32,
                                              const secp256k1_pubkey *pubkeys,
                                              size_t count) {
    int ok;
    secp256k1_ecdsa_recoverable_batch_workspace *workspace;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(signatures != NULL);
    ARG_CHECK(messages32 != NULL);
    ARG_CHECK(pubkeys != NULL);
    if (count == 0) return 1;

    workspace = secp256k1_ecdsa_recoverable_batch_workspace_create(ctx, count);
    if (workspace == NULL) return 0;
    ok = secp256k1_ecdsa_recoverable_verify_batch_workspace(
        ctx, workspace, signatures, messages32, pubkeys, count);
    secp256k1_ecdsa_recoverable_batch_workspace_destroy(ctx, workspace);
    return ok;
}

#endif /* SECP256K1_MODULE_RECOVERY_MAIN_H */
