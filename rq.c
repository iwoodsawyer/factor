/*
 * RQ factorization
 *
 * R = rq(A)
 * [R,Q] = rq(A)
 * [R,Q] = rq(A,0)
 *
 * example compile command (see also make_factor.m):
 * mex -O rq.c libmwlapack.lib
 * or
 * mex -O rq.c libmwblas.lib libmwlapack.lib (>= R2007B)
 *
 * calls the SGERQF/DGERQF/CGERQF/ZGERQF and
 * SORGRQ/DORGRQ/CUNGRQ/ZUNGRQ named LAPACK functions
 */

#include "mex.h"
#include "factor.h"
#include "matrix.h"

void rq_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, lda, m2, dc = 1; 
    mwSignedIndex i, j, start, limit;
    size_t element_size = sizeof(double);
    double *Qpr, *Rpr, *Ipr, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    double *Qpi, *Rpi, *Ipi;
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
                Qpr = mxGetData(plhs[0]);
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
        zgerqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    else {
        dgerqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the DGERQF function */
    if (cplx) {
        zgerqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    else {
        dgerqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("ZGERQF not successful");
        }
        else {
            mexErrMsgTxt("DGERQF not successful");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m < n) {
        plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
            for (j=0; j<m; j++) {
                limit = j<m-1 ? j+1:m;
                for (i=0; i<limit; i++) {
                    Rpr[j*m+i] = Ap[(j+n-m)*2*lda+2*i];
                    Rpi[j*m+i] = Ap[(j+n-m)*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<m; j++) {
                limit = j<m-1 ? j+1:m;
                for (i=0; i<dc*limit; i++) {
                    Rpr[j*dc*m+i] = Ap[(j+n-m)*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
            if (m < n) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<limit; i++) {
                        Rpr[(j+n-m)*m+i] = Ap[(j+n-m)*2*lda+2*i];
                        Rpi[(j+n-m)*m+i] = Ap[(j+n-m)*2*lda+2*i+1];
                    }
                }
            }
            else {
                for (j=0; j<n; j++) {
                    limit = j < min_mn-1 ? j+1 : min_mn;
                    for (i=0; i<limit+m-n; i++) {
                        Rpr[j*m+i] = Ap[j*2*lda+2*i];
                        Rpi[j*m+i] = Ap[j*2*lda+2*i+1];
                    }
                }
            }
        }
        else {
        #endif
            if (m < n) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<dc*limit; i++) {
                        Rpr[(j+n-m)*dc*m+i] = Ap[(j+n-m)*dc*lda+i];
                    }
                }
            }
            else {
                for (j=0; j<n; j++) {
                    limit = j < min_mn-1 ? j+1 : min_mn;
                    for (i=0; i<dc*(limit+m-n); i++) {
                        Rpr[j*dc*m+i] = Ap[j*dc*lda+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX) 
        }
        #endif
    }

    if (nlhs == 2) {
        m2 = (econ == 1) ? min_mn : n;
        
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            if (m > n) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=start+m-n; i<m; i++) {
                        Ap[j*2*lda+2*i-2*m+2*n] = Ap[j*2*lda+2*i];
                        Ap[j*2*lda+2*i-2*m+2*n+1] = Ap[j*2*lda+2*i+1];
                    }
                }
            }
            else if ((m < n) && (econ != 1)) {
                for (j=0; j<n-1; j++) {
                    start = j<n-m ? 0:j-n+m;
                    for (i=m; i>start; i--) {
                        Ap[j*2*lda+2*(i-1)+2*n-2*m] = Ap[j*2*lda+2*(i-1)];
                        Ap[j*2*lda+2*(i-1)+2*n-2*m+1] = Ap[j*2*lda+2*(i-1)+1];
                    }
                }
            }
        }
        else {
        #endif
            if (m > n) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=dc*(start+m-n); i<m; i++) {
                        Ap[j*dc*lda+dc*(n-m)+i] = Ap[j*dc*lda+i];
                    }
                }
            }
            else if ((m < n) && (econ != 1)) {
                for (j=0; j<n-1; j++) {
                    start = j<n-m ? 0:j-n+m;
                    for (i=dc*m; i>dc*start; i--) {
                        Ap[j*dc*lda+dc*(n-m-1)+i] = Ap[j*dc*lda-dc*1+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        
        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            zungrq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        else {
            dorgrq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the DORGRQ function */
        if (cplx) {
            zungrq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        else {
            dorgrq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("ZUNGRQ not successful");
            }
            else {
                mexErrMsgTxt("DORGRQ not successful");
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


void rq_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, lda, m2, dc = 1; 
    mwSignedIndex i, j, start, limit;
    size_t element_size = sizeof(float);
    float *Qpr, *Rpr, *Ipr, *Ap, *ptau, *pwork, *pwork2, *psize, *psize2;
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    float *Qpi, *Rpi, *Ipi;
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
                Qpr = mxGetData(plhs[0]);
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
        cgerqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    else {
        sgerqf(&m, &n, Ap, &lda, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the DGERQF function */
    if (cplx) {
        cgerqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    else {
        sgerqf(&m, &n, Ap, &lda, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("CGERQF not successful");
        }
        else {
            mexErrMsgTxt("SGERQF not successful");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m < n) {
        plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
            for (j=0; j<m; j++) {
                limit = j<m-1 ? j+1:m;
                for (i=0; i<limit; i++) {
                    Rpr[j*m+i] = Ap[(j+n-m)*2*lda+2*i];
                    Rpi[j*m+i] = Ap[(j+n-m)*2*lda+2*i+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<m; j++) {
                limit = j<m-1 ? j+1:m;
                for (i=0; i<dc*limit; i++) {
                    Rpr[j*dc*m+i] = Ap[(j+n-m)*dc*lda+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
        Rpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Rpi = mxGetImagData(plhs[0]);
            if (m < n) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<limit; i++) {
                        Rpr[(j+n-m)*m+i] = Ap[(j+n-m)*2*lda+2*i];
                        Rpi[(j+n-m)*m+i] = Ap[(j+n-m)*2*lda+2*i+1];
                    }
                }
            }
            else {
                for (j=0; j<n; j++) {
                    limit = j < min_mn-1 ? j+1 : min_mn;
                    for (i=0; i<limit+m-n; i++) {
                        Rpr[j*m+i] = Ap[j*2*lda+2*i];
                        Rpi[j*m+i] = Ap[j*2*lda+2*i+1];
                    }
                }
            }
        }
        else {
        #endif
            if (m < n) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<dc*limit; i++) {
                        Rpr[(j+n-m)*dc*m+i] = Ap[(j+n-m)*dc*lda+i];
                    }
                }
            }
            else {
                for (j=0; j<n; j++) {
                    limit = j < min_mn-1 ? j+1 : min_mn;
                    for (i=0; i<dc*(limit+m-n); i++) {
                        Rpr[j*dc*m+i] = Ap[j*dc*lda+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    if (nlhs == 2) {
        m2 = (econ == 1) ? min_mn : n;
        
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            if (m > n) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=start+m-n; i<m; i++) {
                        Ap[j*2*lda+2*i-2*m+2*n] = Ap[j*2*lda+2*i];
                        Ap[j*2*lda+2*i-2*m+2*n+1] = Ap[j*2*lda+2*i+1];
                    }
                }
            }
            else if ((m < n) && (econ != 1)) {
                for (j=0; j<n-1; j++) {
                    start = j<n-m ? 0:j-n+m;
                    for (i=m; i>start; i--) {
                        Ap[j*2*lda+2*(i-1)+2*n-2*m] = Ap[j*2*lda+2*(i-1)];
                        Ap[j*2*lda+2*(i-1)+2*n-2*m+1] = Ap[j*2*lda+2*(i-1)+1];
                    }
                }
            }
        }
        else {
        #endif
            if (m > n) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=dc*(start+m-n); i<dc*m; i++) {
                        Ap[j*dc*lda+dc*(n-m)+i] = Ap[j*dc*lda+i];
                    }
                }
            }
            else if ((m < n) && (econ != 1)) {
                for (j=0; j<n-1; j++) {
                    start = j<n-m ? 0:j-n+m;
                    for (i=dc*m; i>dc*start; i--) {
                        Ap[j*dc*lda+dc*(n-m-1)+i] = Ap[j*dc*lda-dc*1+i];
                    }
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
        
        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            cungrq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        else {
            sorgrq(&m2, &n, &min_mn, Ap, &lda, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the DORGRQ function */
        if (cplx) {
            cungrq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        else {
            sorgrq(&m2, &n, &min_mn, Ap, &lda, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("CUNGRQ not successful");
            }
            else {
                mexErrMsgTxt("SORGRQ not successful");
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
        mexErrMsgTxt("RQ requires one or two input arguments.");
    }
    if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (mxIsDouble(prhs[0])) {
        rq_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0])) {
        rq_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "Class is not supported." );
    }
}
