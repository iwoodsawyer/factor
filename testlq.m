n = 211;
m = 211;

A = randn(m,n);

tic
[L,Q] = lq(A);
toc
norm(A-L*Q)

tic
[L,Q] = lq(A,0);
toc
norm(A-L*Q)

n = 300;
m = 411;

A = randn(m,n);

tic
[L,Q] = lq(A);
toc
norm(A-L*Q)

tic
[L,Q] = lq(A,0);
toc
norm(A-L*Q)

n = 422;
m = 311;

A = randn(m,n);

tic
[L,Q] = lq(A);
toc
norm(A-L*Q)

tic
[L,Q] = lq(A,0);
toc
norm(A-L*Q)