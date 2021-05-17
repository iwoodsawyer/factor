[![View QR/RQ/QL/LQ factorizations on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/16536-qr-rq-ql-lq-factorizations)

Currently RQ, QL, and LQ factorizations are not included in Matlab, although these factorizations can also be done by QR function and additional matrix manipulations if matrix is square. Therefore I wrote these mex files, which uses the internal LAPACK routines of Matlab. QR1 is added to complete the set. They can also handle empty matrices. Enforcing positive elements on diagonal R matrix or column pivoting is supported by the QR1 factorization only.


