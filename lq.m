%LQ     Orthogonal-triangular decomposition.
%   [L,Q] = LQ(A), where A is m-by-n, produces an m-by-n lower triangular
%   matrix L and an n-by-n unitary matrix Q so that A = L*Q.
%
%   [L,Q] = LQ(A,0) produces the "economy size" decomposition.
%   If n>m, only the first m rows of Q and the first m columns of L are
%   computed. If n<=m, this is the same as [L,Q] = LQ(A).
%
%   X = LQ(A) and X = LQ(A,0) return the output of LAPACK's *GELQF
%   routine. TRIU(X) is the lower triangular factor L.
%
%   See also QR.