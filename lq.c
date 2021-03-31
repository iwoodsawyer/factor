/*
 * LQ factorization
 *
 * L = lq(A)
 * [L,Q] = lq(A)
 * [L,Q] = lq(A,0)
 *
 * example compile command (see also make_factor.m):
 * mex -O lq.c libmwlapack.lib
 * or
 * mex -O lq.c libmwblas.lib libmwlapack.lib (>= R2007B)
 *
 * calls the SGELQF/DGELQF/CGELQF/ZGELQF and
 * SORGLQ/DORGLQ/CUNGLQ/ZUNGLQ named LAPACK functions
 */

#include "mex.h"
#include "factor.h"
#include "matrix.h"

void lq_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, lda, m2, dc = 1;
    mwSignedIndex i, j, start;
    size_t element_size = sizeof(double);
    double *Qpr, *Lpr, *Ipr, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    double *Qpi, *Lpi, *Ipi;
    #endif
    mxClassID classid = mxDOUBLE_CLASS;
    mxComplexity cplxflag = mxREAL;

    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }

    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (n != 0) {
                Qpr = mxGetData(plhs[1]);
                for (j=0; j<n; j++) {
                    Qpr[j*dc*n + dc*j] = 1;
                }
            }
        }
        return;
    }
    if (nrhs == 2) {
        if (mxGetScalar(prhs[1]) == 0) {
            econ = 1;
        }
    }

    /* allocate tau */
    min_mn = min(m,n);
    ptau = mxMalloc(dc*min_mn*element_size);

    /* allocate A matrix */
    if ((m < n) && (econ != 1)) {
        lda = n;
    }
    else {
        lda = m;
    }
    Ap = mxMalloc(dc*lda*n*element_size);
    
    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*lda+2*i] = Ipr[j*m+i];
                Ap[j*2*lda+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        if (lda == m) {
            memcpy(Ap,Ipr,dc*m*n*element_size);
        }
        else {
            for (j=0; j<n; j++) {
                for (i=0; i<dc*m; i++) {
                    Ap[j*dc*lda+i] = Ipr[j*dc*m+i];
                }
            }
        }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif

    /* determine blocksize */
    lwork = -1;
    psize = mxMalloc(dc*element_size);
    if (cplx) {
        zgelqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    else {
        dgelqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the DGELQF function */
    if (cplx) {
        zgelqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    else {
        dgelqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("ZGELQF not successful");
        }
        else {
            mexErrMsgTxt("DGELQF not successful");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m < n) {
        plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
        Lpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Lpi = mxGetImagData(plhs[0]);
            for (j=0; j<m; j++) {
                start = j<m ? j:m;
                for (i=start; i<m; i++) {
                    Lpr[j*m+i] = Ap[j*2*lda+2*i];
                    Lpi[j*m+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<m; j++) {
                start = j<m ? j:m;
                for (i=dc*start; i<dc*m; i++) {
                    Lpr[j*dc*m+i] = Ap[j*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        Lpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Lpi = mxGetImagData(plhs[0]);
            for (j=0; j<n; j++) {
                start = j<min_mn ? j:min_mn;
                for (i=start; i<m; i++) {
                    Lpr[j*m+i] = Ap[j*2*lda+2*i];
                    Lpi[j*m+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<n; j++) {
                start = j<min_mn ? j:min_mn;
                for (i=dc*start; i<dc*m; i++) {
                    Lpr[j*dc*m+i] = Ap[j*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    if (nlhs == 2) {
        m2 = (econ == 1) ? min_mn : n;

        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            zunglq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        else {
            dorglq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the DORGLQ function */
        if (cplx) {
            zunglq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        else {
            dorglq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("ZUNGLQ not successful");
            }
            else {
                mexErrMsgTxt("DORGLQ not successful");
            }
        }

        /* allocate Q matrix */
        plhs[1] = mxCreateNumericMatrix(m2,n,classid,cplxflag);
        Qpr = mxGetData(plhs[1]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Qpi = mxGetImagData(plhs[1]);

            /* copy Ap to Qp */
            for (j=0; j<n; j++) {
                for (i=0; i<m2; i++) {
                    Qpr[j*m2+i] = Ap[j*2*lda+2*i];
                    Qpi[j*m2+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            /* copy Ap to Qp */
            if (lda == m2) {
                memcpy(Qpr,Ap,dc*m2*n*element_size);
            }
            else {
                for (j=0; j<n; j++) {
                    for (i=0; i<dc*m2; i++) {
                        Qpr[j*dc*m2+i] = Ap[j*dc*lda+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    mxFree(ptau);
    mxFree(Ap);
}


void lq_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, lda, m2, dc = 1;
    mwSignedIndex i, j, start;
    size_t element_size = sizeof(float);
    float *Qpr, *Lpr, *Ipr, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    float *Qpi, *Lpi, *Ipi;
    #endif
    mxClassID classid = mxSINGLE_CLASS;
    mxComplexity cplxflag = mxREAL;

    /* check complex */
    if (mxIsComplex(prhs[0])) {
        cplxflag = mxCOMPLEX;
        cplx = 1;
        dc = 2;
    }

    /* get matrix data */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (m == 0 || n == 0) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        }
        else {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            if (n != 0) {
                Qpr = mxGetData(plhs[1]);
                for (j=0; j<n; j++) {
                    Qpr[j*dc*n + dc*j] = 1;
                }
            }
        }
        return;
    }
    if (nrhs == 2) {
        if (mxGetScalar(prhs[1]) == 0) {
            econ = 1;
        }
    }

    /* allocate tau */
    min_mn = min(m,n);
    ptau = mxMalloc(dc*min_mn*element_size);

    /* allocate A matrix */
    if ((m < n) && (econ != 1)) {
        lda = n;
    }
    else {
        lda = m;
    }
    Ap = mxMalloc(dc*lda*n*element_size);
    
    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*lda+2*i] = Ipr[j*m+i];
                Ap[j*2*lda+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        if (lda == m) {
            memcpy(Ap,Ipr,dc*m*n*element_size);
        }
        else {
            
            for (j=0; j<n; j++) {
                for (i=0; i<dc*m; i++) {
                    Ap[j*dc*lda+i] = Ipr[j*dc*m+i];
                }
            }
        }
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif

    /* determine blocksize */
    lwork = -1;
    psize = mxMalloc(dc*element_size);
    if (cplx) {
        cgelqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    else {
        sgelqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the SGELQF function */
    if (cplx) {
        cgelqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    else {
        sgelqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("CGELQF not successful");
        }
        else {
            mexErrMsgTxt("SGELQF not successful");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m < n) {
        plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
        Lpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Lpi = mxGetImagData(plhs[0]);
            for (j=0; j<m; j++) {
                start = j<m ? j:m;
                for (i=start; i<m; i++) {
                    Lpr[j*m+i] = Ap[j*2*lda+2*i];
                    Lpi[j*m+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<m; j++) {
                start = j<m ? j:m;
                for (i=start; i<dc*m; i++) {
                    Lpr[j*dc*m+i] = Ap[j*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        Lpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Lpi = mxGetImagData(plhs[0]);
            for (j=0; j<n; j++) {
                start = j<min_mn ? j:min_mn;
                for (i=start; i<m; i++) {
                    Lpr[j*m+i] = Ap[j*2*lda+2*i];
                    Lpi[j*m+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<n; j++) {
                start = j<min_mn ? j:min_mn;
                for (i=dc*start; i<dc*m; i++) {
                    Lpr[j*dc*m+i] = Ap[j*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    if (nlhs == 2) {
        m2 = (econ == 1) ? min_mn : n;

        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            cunglq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        else {
            sorglq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the SORGLQ function */
        if (cplx) {
            cunglq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        else {
            sorglq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("CUNGLQ not successful");
            }
            else {
                mexErrMsgTxt("SORGLQ not successful");
            }
        }

        /* allocate Q matrix */
        plhs[1] = mxCreateNumericMatrix(m2,n,classid,cplxflag);
        Qpr = mxGetData(plhs[1]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Qpi = mxGetImagData(plhs[1]);

            /* copy Ap to Qp */
            for (j=0; j<n; j++) {
                for (i=0; i<m2; i++) {
                    Qpr[j*m2+i] = Ap[j*2*lda+2*i];
                    Qpi[j*m2+i] = Ap[j*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            /* copy Ap to Qp */
            if (lda == m2) {
                memcpy(Qpr,Ap,dc*m2*n*element_size);
            }
            else {
                for (j=0; j<n; j++) {
                    for (i=0; i<dc*m2; i++) {
                        
                        Qpr[j*dc*m2+i] = Ap[j*dc*lda+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    mxFree(ptau);
    mxFree(Ap);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* check for proper number of arguments */
    if (nrhs != 1 && nrhs != 2) {
        mexErrMsgTxt("LQ requires one or two input arguments.");
    }
    if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (mxIsDouble(prhs[0])) {
        lq_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0])) {
        lq_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "Class is not supported." );
    }
}

