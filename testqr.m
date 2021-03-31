n = 211;
m = 211;

A = randn(m,n);

tic
[Q,R] = qr(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr(A,0);
toc
norm(A-Q*R)

tic
[Q,R,E] = qr(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr(A,'vector');
toc
norm(A(:,e)-Q*R)


n = 300;
m = 411;

A = randn(m,n);

tic
[Q,R] = qr(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr(A,0);
toc
norm(A-Q*R)

tic
[Q,R,E] = qr(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr1(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr(A,'vector');
toc
norm(A(:,e)-Q*R)



n = 422;
m = 311;

A = randn(m,n);

tic
[Q,R] = qr(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr(A,0);
toc
norm(A-Q*R)

tic
[Q,R,E] = qr(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr(A,'vector');
toc
norm(A(:,e)-Q*R)