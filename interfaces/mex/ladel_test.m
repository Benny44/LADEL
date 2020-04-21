%% Demo on how to use LADEL

n = 10;
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2;
x = rand(n,1);

%% Example 1: factorize and solve
solver = ladel(n);
solver.factorize(M);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

%% Example 2: factorize_advanced and rowmod
Mbasis = sprand(n,n, 8e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;
Mbasis = Mbasis + M; %make sure entries of M are in Mbasis

solver.factorize_advanced(M, Mbasis);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.delete();