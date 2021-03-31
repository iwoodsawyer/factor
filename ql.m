%QL     Orthogonal-triangular decomposition.
%   [Q,L] = QL(A), where A is m-by-n, produces an m-by-n lower triangular
%   matrix L and an m-by-m unitary matrix Q so that A = Q*L.
%
%   [Q,L] = QL(A,0) produces the "economy size" decomposition.
%   If m>n, only the last n columns of Q and the last n rows of L are
%   computed. If m<=n, this is the same as [Q,L] = QL(A).
%
%   X = QL(A) and X = QL(A,0) return the output of LAPACK's *GEQLF
%   routine. TRIU(X) is the lower triangular factor L.
%
%   See also QR.