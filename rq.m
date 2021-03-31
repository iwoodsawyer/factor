%RQ     Orthogonal-triangular decomposition.
%   [R,Q] = RQ(A), where A is m-by-n, produces an m-by-n upper triangular
%   matrix R and an n-by-n unitary matrix Q so that A = R*Q.
%
%   [R,Q] = RQ(A,0) produces the "economy size" decomposition.
%   If n>m, only the last m rows of Q and the last m columns of R are
%   computed. If n<=m, this is the same as [R,Q] = RQ(A).
%
%   X = RQ(A) and X = RQ(A,0) return the output of LAPACK's *GERQF
%   routine. TRIU(X) is the upper triangular factor R.
%
%   See also QR.