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

static int secp256k1_ecdsa_batch_callback(secp256k1_scalar *scalar, secp256k1_ge *point,
                                           size_t idx, void *data) {
    secp256k1_ecdsa_batch_data *batch = (secp256k1_ecdsa_batch_data *)data;
    *scalar = batch->scalars[idx];
    *point = batch->points[idx];
    return 1;
}

int secp256k1_ecdsa_recoverable_verify_batch(const secp256k1_context *ctx,
                                              const secp256k1_ecdsa_recoverable_signature *signatures,
                                              const unsigned char *messages32,
                                              const secp256k1_pubkey *pubkeys,
                                              size_t count) {
    static const unsigned char domain[] = "Harbor/K1BatchVerify/v1";
    secp256k1_ecdsa_batch_data batch;
    secp256k1_scalar g_scalar = secp256k1_scalar_zero;
    secp256k1_gej result;
    secp256k1_scratch *scratch = NULL;
    secp256k1_sha256 transcript;
    unsigned char transcript_hash[32];
    size_t i;
    int ok = 0;

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(signatures != NULL);
    ARG_CHECK(messages32 != NULL);
    ARG_CHECK(pubkeys != NULL);
    if (count == 0) return 1;
    if (count > SIZE_MAX / (2 * sizeof(*batch.scalars))) return 0;

    batch.scalars = (secp256k1_scalar *)malloc(2 * count * sizeof(*batch.scalars));
    batch.points = (secp256k1_ge *)malloc(2 * count * sizeof(*batch.points));
    if (batch.scalars == NULL || batch.points == NULL) goto cleanup;

    secp256k1_sha256_initialize(&transcript);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, domain, sizeof(domain) - 1);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, (const unsigned char *)signatures, count * sizeof(*signatures));
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, messages32, count * 32);
    secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &transcript, (const unsigned char *)pubkeys, count * sizeof(*pubkeys));
    secp256k1_sha256_finalize(secp256k1_get_hash_context(ctx), &transcript, transcript_hash);

    for (i = 0; i < count; ++i) {
        secp256k1_scalar r, s, z, coefficient, term;
        secp256k1_fe x;
        secp256k1_sha256 coefficient_hash;
        unsigned char coefficient_bytes[32];
        unsigned char r_bytes[32];
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

        secp256k1_scalar_get_b32(r_bytes, &r);
        if (!secp256k1_fe_set_b32_limit(&x, r_bytes)) goto cleanup;
        if (recid & 2) {
            if (secp256k1_fe_cmp_var(&x, &secp256k1_ecdsa_const_p_minus_order) >= 0) goto cleanup;
            secp256k1_fe_add(&x, &secp256k1_ecdsa_const_order_as_fe);
        }
        if (!secp256k1_ge_set_xo_var(&batch.points[2 * i], &x, recid & 1)) goto cleanup;
        if (!secp256k1_pubkey_load(ctx, &batch.points[2 * i + 1], &pubkeys[i])) goto cleanup;

        for (j = 0; j < sizeof(index_bytes); ++j)
            index_bytes[sizeof(index_bytes) - 1 - j] = (unsigned char)(i >> (8 * j));
        secp256k1_sha256_initialize(&coefficient_hash);
        secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_hash, transcript_hash, sizeof(transcript_hash));
        secp256k1_sha256_write(secp256k1_get_hash_context(ctx), &coefficient_hash, index_bytes, sizeof(index_bytes));
        secp256k1_sha256_finalize(secp256k1_get_hash_context(ctx), &coefficient_hash, coefficient_bytes);
        secp256k1_scalar_set_b32(&coefficient, coefficient_bytes, &overflow);
        if (overflow || secp256k1_scalar_is_zero(&coefficient))
            secp256k1_scalar_set_int(&coefficient, 1);

        secp256k1_scalar_mul(&batch.scalars[2 * i], &coefficient, &s);
        secp256k1_scalar_mul(&batch.scalars[2 * i + 1], &coefficient, &r);
        secp256k1_scalar_negate(&batch.scalars[2 * i + 1], &batch.scalars[2 * i + 1]);
        secp256k1_scalar_set_b32(&z, messages32 + i * 32, NULL);
        secp256k1_scalar_mul(&term, &coefficient, &z);
        secp256k1_scalar_add(&g_scalar, &g_scalar, &term);
    }
    secp256k1_scalar_negate(&g_scalar, &g_scalar);

    scratch = secp256k1_scratch_create(&ctx->error_callback, 1024 * 1024 + count * 4096);
    if (scratch == NULL) goto cleanup;
    if (!secp256k1_ecmult_multi_var(&ctx->error_callback, scratch, &result, &g_scalar,
                                    secp256k1_ecdsa_batch_callback, &batch, 2 * count)) goto cleanup;
    ok = secp256k1_gej_is_infinity(&result);

cleanup:
    if (scratch != NULL) secp256k1_scratch_destroy(&ctx->error_callback, scratch);
    free(batch.scalars);
    free(batch.points);
    return ok;
}

#endif /* SECP256K1_MODULE_RECOVERY_MAIN_H */
