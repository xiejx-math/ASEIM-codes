% This Matlab file is used to compare REABK, AREABK, and AmREABK
% This file can reproduce Table 1 in the manuscript

close all;
clear;

%%
run_times=20; % average times

%% generated the matrix A using the data from SuiteSparse Matrix Collection

%load ash958;
load WorldCities;
%load Franz1.mat;
%load crew1.mat;
%load model1.mat;
%load bibd_16_8.mat; 
A=Problem.A;
[m,n]=size(A);

ell=30; % block size

%% some vectors are used to store the desired numerical results
CPU_REABK=zeros(run_times,1);
CPU_AREABK=zeros(run_times,1);
CPU_AmREABK=zeros(run_times,1);

Iter_REABK=zeros(run_times,1);
Iter_AREABK=zeros(run_times,1);
Iter_AmREABK=zeros(run_times,1);

for jj=1:run_times

    %% generated the right-hand vector b

    x=randn(n,1);
    Z=null(full(A)');
    r=Z*randn(size(Z,2),1);
    b=A*x+r;
    xLS=lsqminnorm(A,b);

    %% parameter setup
    opts.xstar=xLS;

    if ell==1
        TOL1=10^(-24); 
        % if the block size is 1, we need to modify the rule to verify whether a varible is equal to zero
    end

    %% AREABK

    [xAREABK,OutAREABK]=My_AREABK(A,b,ell,opts);

    %% AmREABK
    [xAmREABK,OutAmREABK]=My_AmREABK(A,b,ell,opts);

    %% REABK
    [xREABK,OutREABK]=My_REABK(A,b,ell,opts);
 
    %% update and store the results

    CPU_REABK(jj)=OutREABK.times(end);
    Iter_REABK(jj)=OutREABK.iter;

    CPU_AREABK(jj)=OutAREABK.times(end);
    Iter_AREABK(jj)=OutAREABK.iter;

    CPU_AmREABK(jj)=OutAmREABK.times(end);
    Iter_AmREABK(jj)=OutAmREABK.iter;
    %%

    if OutAmREABK.iter>=2000000
        OutAmREABK.iter
        OutAmREABK.times(end)
        fprintf('Max_iter\n');
        return;
    end

    %%%%%%
    fprintf('Number of iterations: %d,%d,%d\n',OutREABK.iter,OutAREABK.iter,OutAmREABK.iter)

end

fprintf('m=%d, n=%d\n',m,n)
fprintf('Average number of iterations: REABK=%8.2f, AREABK=%8.2f, AmREABK=%8.2f\n',mean(Iter_REABK),mean(Iter_AREABK),mean(Iter_AmREABK))

fprintf('Average CPU time: REABK=%8.4f, AREABK=%8.4f, AmREABK=%8.4f\n',mean(CPU_REABK),mean(CPU_AREABK),mean(CPU_AmREABK))

fprintf('%8.2f & %8.4f & %8.2f & %8.4f & %8.2f & %8.4f\n',mean(Iter_REABK),mean(CPU_REABK),...
    mean(Iter_AREABK),mean(CPU_AREABK),mean(Iter_AmREABK),mean(CPU_AmREABK))



