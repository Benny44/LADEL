%% Demo on how to use LADEL

n = 100;
M = sprand(n,n, 1e-1);

solver = ladel(n);




solver.delete();