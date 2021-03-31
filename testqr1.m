n = 211;
m = 211;

A = randn(m,n);

tic
[Q,R] = qr1(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,'pos');
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0,'pos');
toc
norm(A-Q*R)

tic
[Q,R,E] = qr1(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr1(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr1(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr1(A,'vector');
toc
norm(A(:,e)-Q*R)


n = 300;
m = 411;

A = randn(m,n);

tic
[Q,R] = qr1(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,'pos');
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0,'pos');
toc
norm(A-Q*R)

tic
[Q,R,E] = qr1(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr1(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr1(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr1(A,'vector');
toc
norm(A(:,e)-Q*R)



n = 422;
m = 311;

A = randn(m,n);

tic
[Q,R] = qr1(A);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0);
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,'pos');
toc
norm(A-Q*R)

tic
[Q,R] = qr1(A,0,'pos');
toc
norm(A-Q*R)

tic
[Q,R,E] = qr1(A);
toc
norm(A*E-Q*R)

tic
[Q,R,E] = qr1(A,'matrix');
toc
norm(A*E-Q*R)

tic
[Q,R,e] = qr1(A,0);
toc
norm(A(:,e)-Q*R)

tic
[Q,R,e] = qr1(A,'vector');
toc
norm(A(:,e)-Q*R)