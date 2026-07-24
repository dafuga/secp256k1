/***********************************************************************
 * Copyright (c) 2026 Harbor contributors                              *
 * Distributed under the MIT software license.                         *
 ***********************************************************************/

#ifndef SECP256K1_HARBOR_BATCH_BACKEND_PORTABLE_IMPL_H
#define SECP256K1_HARBOR_BATCH_BACKEND_PORTABLE_IMPL_H

/* Portable reference backend. Independent addition chains execute in
 * lockstep, which also gives architecture backends a differential oracle. */
static int secp256k1_harbor_batch_fe_sqrt(
        secp256k1_fe r[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH],
        const secp256k1_fe a[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH]) {
    secp256k1_fe x2[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x3[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x6[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x9[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x11[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x22[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x44[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x88[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x176[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x220[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe x223[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    secp256k1_fe t1[HARBOR_SECP256K1_BATCH_BACKEND_WIDTH];
    size_t i;
    int j;
    int ret = 1;

    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        VERIFY_CHECK(&r[i] != &a[i]);
        SECP256K1_FE_VERIFY(&a[i]);
        SECP256K1_FE_VERIFY_MAGNITUDE(&a[i], 8);
        secp256k1_fe_sqr(&x2[i], &a[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x2[i], &x2[i], &a[i]);
        secp256k1_fe_sqr(&x3[i], &x2[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x3[i], &x3[i], &a[i]);
        x6[i] = x3[i];
    }
    for (j = 0; j < 3; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x6[i], &x6[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x6[i], &x6[i], &x3[i]);
        x9[i] = x6[i];
    }
    for (j = 0; j < 3; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x9[i], &x9[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x9[i], &x9[i], &x3[i]);
        x11[i] = x9[i];
    }
    for (j = 0; j < 2; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x11[i], &x11[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x11[i], &x11[i], &x2[i]);
        x22[i] = x11[i];
    }
    for (j = 0; j < 11; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x22[i], &x22[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x22[i], &x22[i], &x11[i]);
        x44[i] = x22[i];
    }
    for (j = 0; j < 22; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x44[i], &x44[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x44[i], &x44[i], &x22[i]);
        x88[i] = x44[i];
    }
    for (j = 0; j < 44; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x88[i], &x88[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x88[i], &x88[i], &x44[i]);
        x176[i] = x88[i];
    }
    for (j = 0; j < 88; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x176[i], &x176[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x176[i], &x176[i], &x88[i]);
        x220[i] = x176[i];
    }
    for (j = 0; j < 44; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x220[i], &x220[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x220[i], &x220[i], &x44[i]);
        x223[i] = x220[i];
    }
    for (j = 0; j < 3; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&x223[i], &x223[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&x223[i], &x223[i], &x3[i]);
        t1[i] = x223[i];
    }
    for (j = 0; j < 23; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&t1[i], &t1[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_mul(&t1[i], &t1[i], &x22[i]);
    for (j = 0; j < 6; ++j) {
        for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) secp256k1_fe_sqr(&t1[i], &t1[i]);
    }
    for (i = 0; i < HARBOR_SECP256K1_BATCH_BACKEND_WIDTH; ++i) {
        secp256k1_fe_mul(&t1[i], &t1[i], &x2[i]);
        secp256k1_fe_sqr(&t1[i], &t1[i]);
        secp256k1_fe_sqr(&r[i], &t1[i]);
        secp256k1_fe_sqr(&t1[i], &r[i]);
        ret &= secp256k1_fe_equal(&t1[i], &a[i]);
    }
    return ret;
}

#endif /* SECP256K1_HARBOR_BATCH_BACKEND_PORTABLE_IMPL_H */
