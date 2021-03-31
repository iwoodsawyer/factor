%QR1     Orthogonal-triangular decomposition.
%   [Q,R] = QR1(A), where A is m-by-n, produces an m-by-n upper triangular
%   matrix R and an m-by-m unitary matrix Q so that A = Q*R.
%
%   [Q,R] = QR1(A,'pos') enforces positive elements on diagonal R matrix. 
%
%   [Q,R] = QR1(A,0) produces the "economy size" decomposition.
%   If m>n, only the first n columns of Q and the first n rows of R are
%   computed. If m<=n, this is the same as [Q,R] = QR1(A).
%
%   [Q,R] = QR1(A,0,'pos') enforces positive elements on diagonal R matrix
%   and produces the "economy size" decomposition. 
%
%   [Q,R,E] = QR1(A) or [Q,R,E] = QR1(A,'matrix') produces unitary Q, upper
%   triangular R and a permutation matrix E so that A*E = Q*R. The column
%   permutation E is chosen so that ABS(DIAG(R)) is decreasing.
%
%   [Q,R,e] = QR1(A,'vector') returns the permutation information as a
%   vector instead of a matrix. That is, e is a row vector such that A(:,e)
%   = Q*R.
%
%   [Q,R,e] = QR1(A,0) produces an economy-size decomposition in which e is
%   a permutation vector, so that A(:,e) = Q*R.
%
%   X = QR1(A) and X = QR1(A,0) return a matrix X such that TRIU(X) is the
%   upper triangular factor R.
%
%   See also QR.