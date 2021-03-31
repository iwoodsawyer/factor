n = 211;
m = 211;

A = randn(m,n);

tic
[Q,L] = ql(A);
toc
norm(A-Q*L)

tic
[Q,L] = ql(A,0);
toc
norm(A-Q*L)

n = 300;
m = 411;

A = randn(m,n);

tic
[Q,L] = ql(A);
toc
norm(A-Q*L)

tic
[Q,L] = ql(A,0);
toc
norm(A-Q*L)

n = 422;
m = 311;

A = randn(m,n);

tic
[Q,L] = ql(A);
toc
norm(A-Q*L)

tic
[Q,L] = ql(A,0);
toc
norm(A-Q*L)