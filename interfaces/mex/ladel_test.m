%% Demo on how to use LADEL
if exist('solver')
    solver.delete();
end

n = 4;
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2;
x = rand(n,1);

%% Example 1: factorize and solve
solver = ladel(n);
solver.factorize(M);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

%% Example 2: factorize_advanced and row_mod
Mbasis = sprand(n,n, 8e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make the n/2 row/column only contain a diagonal element
M(n/2,:) = zeros(1,n);
M(:,n/2) = zeros(n,1);
d = rand(1);
M(n/2,n/2) = 1; 

Mbasis = Mbasis + M; %make sure entries of M are in Mbasis

solver.factorize_advanced(M, Mbasis);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.row_mod(n/2, Mbasis(:,n/2), full(Mbasis(n/2,n/2)));
Mupd = M;
Mupd(:,n/2) = Mbasis(:,n/2);
Mupd(n/2,:) = Mbasis(n/2,:);

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < 1e-12);




solver.delete();