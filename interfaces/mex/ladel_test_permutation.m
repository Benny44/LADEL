%% Demo on how to use LADEL (without permutation)
if exist('solver')
    solver.delete();
    clear solver
end

test_cholmod = true;

ordering = 1; %1 for AMD, 0 for natural ordering

tic
n = 1000;
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2;
x = rand(n,1);
times.generate.M = toc;

%% Example 1: factorize and solve
solver = ladel(n);
% [L,D,p] = solver.factorize(M, ordering);
solver.factorize(M, ordering);

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

%% Example 2: factorize_advanced and row_mod
tic;
Mbasis = sprand(n,n, 2e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make the n/2 row/column only contain a diagonal element
M(n/2,:) = zeros(1,n);
M(:,n/2) = zeros(n,1);
d = rand(1);
M(n/2,n/2) = 1; 

Mbasis = Mbasis + M; %make sure entries of M are in Mbasis
times.generate.Mbasis = toc;

% clear all
% load('problematic');
% clear times
% clear solver
% solver = ladel(n);
% ordering = 1;
tic;
% [L,D,p] = solver.factorize_advanced(M, Mbasis, ordering);
solver.factorize_advanced(M, Mbasis, ordering);
times.ladel.factorize = toc;
if test_cholmod
    tic;
    [LD,~, LD_p] = ldlchol(M);
    LD_pinv = 1:n;
    LD_pinv(LD_p) = LD_pinv;
    times.chol.factorize = toc;
end

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);
if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < 1e-12);
end

%% ADD row using row_mod
clear row;
row = Mbasis(:,n/2);
tic;
% [Lupd,Dupd,p] = solver.row_mod(n/2, row, full(Mbasis(n/2,n/2)));
solver.row_mod(n/2, row, full(Mbasis(n/2,n/2)));
times.ladel.rowadd = toc;

if test_cholmod
    row = Mbasis(:,n/2);
    tic;
    row = row(LD_p);
    LD = ldlrowmod(LD, LD_pinv(n/2), row);
    times.chol.rowadd = toc;
end

Mupd = M;
Mupd(:,n/2) = Mbasis(:,n/2);
Mupd(n/2,:) = Mbasis(n/2,:);

if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-Mupd\x) < 1e-12);
end

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < 1e-12);

%% DELETE row using row_mod
% We delete the row we just added again

if test_cholmod
    tic;
    LD = ldlrowmod(LD, LD_pinv(n/2));
    times.chol.rowdel = toc;
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < 1e-12);
end

tic;
% [Lfinal, Dfinal, p] = solver.row_mod(n/2);
solver.row_mod(n/2);
times.ladel.rowdel = toc;

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.delete();

times.ladel
times.chol