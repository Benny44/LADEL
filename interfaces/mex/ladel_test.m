%% Demo on how to use LADEL

n = 10;
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2;
% M = rand(n,n);
x = rand(n,1);

solver = ladel(n);
solver.factorize(M);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.delete();