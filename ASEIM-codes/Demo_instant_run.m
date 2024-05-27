% This Matlab file provided is an instant run code for AmREABK, AREABK, and REABK 
close all;
clear;

m=1024;
n=128;

rank=100;

kappa=40;

[U,~]=qr(randn(m, rank), 0);
[V,~]=qr(randn(n, rank), 0);
D = diag(1+(kappa-1).*rand(rank, 1));
A=U*D*V';

%% generated the right-hand vector b

A=full(A);
x=randn(n,1);
Z=null(full(A)');
r=Z*randn(size(Z,2),1);
b=A*x+r;
if rank==n
    xLS=x;
else
    xLS=lsqminnorm(A,b);
end

%% parameter setup
opts.xstar=xLS;
%%
ell=50;
%%
[xAREABK,OutAREABK]=My_AREABK(A,b,ell,opts);

%%
[xAmREABK,OutAmREABK]=My_AmREABK(A,b,ell,opts);

%% REABK
[xREABK,OutREABK]=My_REABK(A,b,ell,opts);

%%
fprintf('m=%d,n=%d,rank=%d,block size=%d, kappa=%d\n',m,n,rank,ell,kappa)
fprintf('AREABK: CPU time %2.4f\n',OutAREABK.times(end))
fprintf('AmREABK: CPU time %2.4f\n',OutAmREABK.times(end))
fprintf('REABK: CPU time %2.4f\n',OutREABK.times(end))

fprintf('AREABK: RSE  %2.4e\n',OutAREABK.error(end))
fprintf('AmREABK: RSE  %2.4e\n',OutAmREABK.error(end))
fprintf('REABK: RSE  %2.4e\n',OutREABK.error(end))

fprintf('AREABK: Iter %2.4f\n',OutAREABK.iter)
fprintf('AmREABK: Iter %2.4f\n',OutAmREABK.iter)
fprintf('REABK: Iter %2.4f\n',OutREABK.iter)







