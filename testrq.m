n = 211;
m = 211;

A = randn(m,n);

tic
[R,Q] = rq(A);
toc
norm(A-R*Q)

tic
[R,Q] = rq(A,0);
toc
norm(A-R*Q)

n = 300;
m = 411;

A = randn(m,n);

tic
[R,Q] = rq(A);
toc
norm(A-R*Q)

tic
[R,Q] = rq(A,0);
toc
norm(A-R*Q)

n = 422;
m = 311;

A = randn(m,n);

tic
[R,Q] = rq(A);
toc
norm(A-R*Q)

tic
[R,Q] = rq(A,0);
toc
norm(A-R*Q)