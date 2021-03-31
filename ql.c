/*
 * QL factorization
 *
 * L = ql(A)
 * [Q,L] = ql(A)
 * [Q,L] = ql(A,0)
 *
 * example compile command (see also make_factor.m):
 * mex -O ql.c libmwlapack.lib
 * or
 * mex -O rq.c libmwblas.lib libmwlapack.lib (>= R2007B)
 *
 * calls the SGEQLF/DGEQLF/CGEQLF/ZGEQLF and
 * SORGQL/DORGQL/CUNGQL/ZUNGQL named LAPACK functions
 */

#include "mex.h"
#include "factor.h"
#include "matrix.h"

void ql_double(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, n2, dc = 1;
    mwSignedIndex i, j, start, limit;
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
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (j=0; j<m; j++) {
                    Qpr[j*dc*m + dc*j] = 1;
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
    if (m > n && (econ != 1)) {
        Ap = mxMalloc(dc*m*m*element_size);
    }
    else {
        Ap = mxMalloc(dc*m*n*element_size);
    }

    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        memcpy(Ap,Ipr,dc*m*n*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif

    /* determine blocksize */
    lwork = -1;
    psize = mxMalloc(dc*element_size);
    if (cplx) {
        zgeqlf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    else {
        dgeqlf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the DGEQLF function */
    if (cplx) {
        zgeqlf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    else {
        dgeqlf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    mxFree(pwork);
    if (info != 0) {
        mxFree(ptau);
        mxFree(Ap);
        if (cplx) {
            mexErrMsgTxt("ZGEQLF not successful");
        }
        else {
            mexErrMsgTxt("DGEQLF not successful");
        }
    }

    /* extract lower triangular part */
    if ((econ == 1) && m > n) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            Lpr = mxGetData(plhs[0]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[0]);
            }
            #endif
        }
        else {
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            Lpr = mxGetData(plhs[1]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[1]);
            }
            #endif
        }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            for (j=0; j<n; j++) {
                start = j<n ? j:n;
                for (i=start; i<n; i++) {
                    Lpr[j*n+i] = Ap[j*2*m+2*(i+m-n)];
                    Lpi[j*n+i] = Ap[j*2*m+2*(i+m-n)+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<n; j++) {
                start = j<n ? j:n;
                for (i=dc*start; i<dc*n; i++) {
                    Lpr[j*dc*n+i] = Ap[j*dc*m+dc*(m-n)+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            Lpr = mxGetData(plhs[0]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[0]);
            }
            #endif
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            Lpr = mxGetData(plhs[1]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[1]);
            }
            #endif
        }
        if (m > n) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=start+m-n; i<m; i++) {
                        Lpr[j*m+i] = Ap[j*2*m+2*i];
                        Lpi[j*m+i] = Ap[j*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=dc*(start+m-n); i<dc*m; i++) {
                        Lpr[j*dc*m+i] = Ap[j*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
        else {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<n; j++) {
                    start = j<(n-m) ? 0:j-(n-m);
                    for (i=start; i<m; i++) {
                        Lpr[j*m+i] = Ap[j*2*m+2*i];
                        Lpi[j*m+i] = Ap[j*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<n; j++) {
                    start = j<(n-m) ? 0:j-(n-m);
                    for (i=dc*start; i<dc*m; i++) {
                        Lpr[j*dc*m+i] = Ap[j*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
    }

    if (nlhs == 2) {
        n2 = (econ == 1) ? min_mn : m;

        if (m < n) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<limit; i++) {
                        Ap[j*2*m+2*i] = Ap[(j+n-m)*2*m+2*i];
                        Ap[j*2*m+2*i+1] = Ap[(j+n-m)*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<dc*limit; i++) {
                        Ap[j*dc*m+i] = Ap[(j+n-m)*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
        else if ((m > n) && (econ != 1)) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=n; j>0; j--) {
                    limit = j-1<n ? m-n+j:m;
                    for (i=0; i<limit; i++) {
                        Ap[(j+m-n-1)*2*m+2*i] = Ap[(j-1)*2*m+2*i];
                        Ap[(j+m-n-1)*2*m+2*i+1] = Ap[(j-1)*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=n; j>0; j--) {
                    limit = j-1<n ? m-n+j:m;
                    for (i=0; i<dc*limit; i++) {
                        Ap[(j+m-n-1)*dc*m+i] = Ap[(j-1)*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }

        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            zungql(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        else {
            dorgql(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the DORGQL function */
        if (cplx) {
            zungql(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        else {
            dorgql(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("ZUNGQL not successful");
            }
            else {
                mexErrMsgTxt("DORGQL not successful");
            }
        }

        /* allocate Q matrix */
        plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
        Qpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Qpi = mxGetImagData(plhs[0]);

            /* copy Ap to Qp */
            for (j=0; j<n2; j++) {
                for (i=0; i<m; i++) {
                    Qpr[j*m+i] = Ap[j*2*m+2*i];
                    Qpi[j*m+i] = Ap[j*2*m+2*i+1];
                }
            }
        }
        else {
        #endif
            /* copy Ap to Qp */
            memcpy(Qpr,Ap,dc*m*n2*element_size);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }

    mxFree(ptau);
    mxFree(Ap);
}


void ql_single(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex lwork, info = 0, econ = 0, cplx = 0;
    mwSignedIndex m, n, min_mn, n2, dc = 1;
    mwSignedIndex i, j, start, limit;
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
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            plhs[0] = mxCreateNumericMatrix(m,m,classid,cplxflag);
            if (m != 0) {
                Qpr = mxGetData(plhs[0]);
                for (j=0; j<m; j++) {
                    Qpr[j*dc*m + dc*j] = 1;
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
    if (m > n && (econ != 1)) {
        Ap = mxMalloc(dc*m*m*element_size);
    }
    else {
        Ap = mxMalloc(dc*m*n*element_size);
    }

    /* copy input to A matrix */
    Ipr = mxGetData(prhs[0]);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    if (cplx) {
        Ipi = mxGetImagData(prhs[0]);
        for (j=0; j<n; j++) {
            for (i=0; i<m; i++) {
                Ap[j*2*m+2*i] = Ipr[j*m+i];
                Ap[j*2*m+2*i+1] = Ipi[j*m+i];
            }
        }
    }
    else {
    #endif
        memcpy(Ap,Ipr,dc*m*n*element_size);
    #if !(MX_HAS_INTERLEAVED_COMPLEX)
    }
    #endif

    /* determine blocksize */
    lwork = -1;
    psize = mxMalloc(dc*element_size);
    if (cplx) {
        cgeqlf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    else {
        sgeqlf(&m, &n, Ap, &m, ptau, psize, &lwork, &info);
    }
    lwork = (mwSignedIndex)psize[0];
    mxFree(psize);

    /* allocate workspace */
    pwork = mxMalloc(dc*lwork*element_size);

    /* calls the SGEQLF function */
    if (cplx) {
        cgeqlf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
    }
    else {
        sgeqlf(&m, &n, Ap, &m, ptau, pwork, &lwork, &info);
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
    if ((econ == 1) && m > n) {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            Lpr = mxGetData(plhs[0]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[0]);
            }
            #endif
        }
        else {
            plhs[1] = mxCreateNumericMatrix(n,n,classid,cplxflag);
            Lpr = mxGetData(plhs[1]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[1]);
            }
            #endif
        }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            for (j=0; j<n; j++) {
                start = j<n ? j:n;
                for (i=start; i<n; i++) {
                    Lpr[j*n+i] = Ap[j*2*m+2*(i+m-n)];
                    Lpi[j*n+i] = Ap[j*2*m+2*(i+m-n)+1];
                }
            }
        }
        else {
        #endif
            for (j=0; j<n; j++) {
                start = j<n ? j:n;
                for (i=dc*start; i<dc*n; i++) {
                    Lpr[j*dc*n+i] = Ap[j*dc*m+dc*(m-n)+i];
                }
            }
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        }
        #endif
    }
    else {
        if (nlhs == 1) {
            plhs[0] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            Lpr = mxGetData(plhs[0]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[0]);
            }
            #endif
        }
        else {
            plhs[1] = mxCreateNumericMatrix(m,n,classid,cplxflag);
            Lpr = mxGetData(plhs[1]);
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                Lpi = mxGetImagData(plhs[1]);
            }
            #endif
        }
        if (m > n) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=start+m-n; i<m; i++) {
                        Lpr[j*m+i] = Ap[j*2*m+2*i];
                        Lpi[j*m+i] = Ap[j*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<n; j++) {
                    start = j<n ? j:n;
                    for (i=dc*(start+m-n); i<dc*m; i++) {
                        Lpr[j*dc*m+i] = Ap[j*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
        else {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<n; j++) {
                    start = j<(n-m) ? 0:j-(n-m);
                    for (i=start; i<m; i++) {
                        Lpr[j*m+i] = Ap[j*2*m+2*i];
                        Lpi[j*m+i] = Ap[j*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<n; j++) {
                    start = j<(n-m) ? 0:j-(n-m);
                    for (i=dc*start; i<dc*m; i++) {
                        Lpr[j*dc*m+i] = Ap[j*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
    }

    if (nlhs == 2) {
        n2 = (econ == 1) ? min_mn : m;

        if (m < n) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<limit; i++) {
                        Ap[j*2*m+2*i] = Ap[(j+n-m)*2*m+2*i];
                        Ap[j*2*m+2*i+1] = Ap[(j+n-m)*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=0; j<m; j++) {
                    limit = j<m-1 ? j+1:m;
                    for (i=0; i<dc*limit; i++) {
                        Ap[j*dc*m+i] = Ap[(j+n-m)*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }
        else if ((m > n) && (econ != 1)) {
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            if (cplx) {
                for (j=n; j>0; j--) {
                    limit = j-1<n ? m-n+j:m;
                    for (i=0; i<limit; i++) {
                        Ap[(j+m-n-1)*2*m+2*i] = Ap[(j-1)*2*m+2*i];
                        Ap[(j+m-n-1)*2*m+2*i+1] = Ap[(j-1)*2*m+2*i+1];
                    }
                }
            }
            else {
            #endif
                for (j=n; j>0; j--) {
                    limit = j-1<n ? m-n+j:m;
                    for (i=0; i<dc*limit; i++) {
                        Ap[(j+m-n-1)*dc*m+i] = Ap[(j-1)*dc*m+i];
                    }
                }
            #if !(MX_HAS_INTERLEAVED_COMPLEX)
            }
            #endif
        }

        /* determine blocksize */
        lwork = -1;
        psize2 = mxMalloc(dc*element_size);
        if (cplx) {
            cungql(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        else {
            sorgql(&m, &n2, &min_mn, Ap, &m, ptau, psize2, &lwork, &info);
        }
        lwork = (mwSignedIndex)psize2[0];
        mxFree(psize2);

        /* allocate workspace */
        pwork2 = mxMalloc(dc*lwork*element_size);

        /* calls the SORGQL function */
        if (cplx) {
            cungql(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        else {
            sorgql(&m, &n2, &min_mn, Ap, &m, ptau, pwork2, &lwork, &info);
        }
        mxFree(pwork2);
        if (info != 0) {
            mxFree(ptau);
            mxFree(Ap);
            if (cplx) {
                mexErrMsgTxt("CUNGQL not successful");
            }
            else {
                mexErrMsgTxt("SORGQL not successful");
            }
        }

        /* allocate Q matrix */
        plhs[0] = mxCreateNumericMatrix(m,n2,classid,cplxflag);
        Qpr = mxGetData(plhs[0]);
        #if !(MX_HAS_INTERLEAVED_COMPLEX)
        if (cplx) {
            Qpi = mxGetImagData(plhs[0]);

            /* copy Ap to Qp */
            for (j=0; j<n2; j++) {
                for (i=0; i<m; i++) {
                    Qpr[j*m+i] = Ap[j*2*m+2*i];
                    Qpi[j*m+i] = Ap[j*2*m+2*i+1];
                }
            }
        }
        else {
        #endif
            /* copy Ap to Qp */
            memcpy(Qpr,Ap,dc*m*n2*element_size);
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
        mexErrMsgTxt("QL requires one or two input arguments.");
    }
    if (nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (!mxIsNumeric(prhs[0]) || mxIsSparse(prhs[0])) {
        mexErrMsgTxt( "Input must be a full matrix." );
    }
    if (mxIsDouble(prhs[0])) {
        ql_double(nlhs, plhs, nrhs, prhs);
    }
    else if (mxIsSingle(prhs[0])) {
        ql_single(nlhs, plhs, nrhs, prhs);
    }
    else {
        mexErrMsgTxt( "Class is not supported." );
    }
}
