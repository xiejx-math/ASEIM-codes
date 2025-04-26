% This Matlab file is used to compare AmREABK, AREABK, and REABK
% based on the data from LIBSVM

clear;
close all;

%%
opts.Max_iter=100; % Number of iterations
ell=300;% size of the block
run_times=20; % average times

%% LIBSVM data
%load a9a; %Max_iter=2500
%A=a9a_inst;
%clear a9a;

%load cod-rna %  Max_iter=30000
%A=cod_rna_inst;
%clear cod_rna;

load ijcnn1 % Max_iter=100
A=ijcnn1_inst;
clear ijcnn1


[m,n]=size(A);

%%
opts.TOL=10^(-32);
opts.TOL1=eps^2;
opts.sparsity=1;
%% the vector is used to store the numerical results
RSE_REABK=zeros(run_times,opts.Max_iter);
RSE_AREABK=zeros(run_times,opts.Max_iter);
RSE_AmREABK=zeros(run_times,opts.Max_iter);

CPU_REABK=zeros(run_times,opts.Max_iter);
CPU_AREABK=zeros(run_times,opts.Max_iter);
CPU_AmREABK=zeros(run_times,opts.Max_iter);

for ii=1:run_times
    x=randn(n,1);
    r=randn(m,1);
    r=r-A*(pinv(full(A))*r);
    b=A*x+r;
    xLS=lsqminnorm(A,b);

    %% parameter setup
    opts.xstar=xLS;

    %% AREABK
    [xAREABK,OutAREABK]=My_AREABK(A,b,ell,opts);

    %% AmREABK
    [xAmREABK,OutAmREABK]=My_AmREABK(A,b,ell,opts);

    %% REABK
    [xREABK,OutREABK]=My_REABK(A,b,ell,opts);

    %%
    RSE_REABK(ii,:)=OutREABK.error(1:opts.Max_iter);%
    RSE_AREABK(ii,:)=OutAREABK.error(1:opts.Max_iter);%
    RSE_AmREABK(ii,:)=OutAmREABK.error(1:opts.Max_iter);


    CPU_REABK(ii,:)=OutREABK.times(1:opts.Max_iter);
    CPU_AREABK(ii,:)=OutAREABK.times(1:opts.Max_iter);
    CPU_AmREABK(ii,:)=OutAmREABK.times(1:opts.Max_iter);

    fprintf('Done,iter=%d\n',ii);
end


%% plot errors
xlabel_i=1:opts.Max_iter;
num_iter_array=xlabel_i';

%% 
lightgray =   [0.8 0.8 0.8];
mediumgray =  [0.6 0.6 0.6];
lightred =    [1 0.9 0.9];
mediumred =   [1 0.6 0.6];
lightgreen =  [0.9 1 0.9];
mediumgreen = [0.6 1 0.6];
lightblue =   [0.9 0.9 1];
mediumblue =  [0.6 0.6 1];
lightmagenta =   [1 0.9 1];
mediummagenta =  [1 0.6 1];

%%
display_names = {'REABK','AREABK','AmREABK'};
arrsIter = {RSE_REABK',RSE_AREABK',RSE_AmREABK'};
num_methods = length(arrsIter);
line_colors = {'red','blue','green'};
minmax_colors = { lightred,lightblue,lightgreen};
quant_colors = { mediumred,mediumblue,mediumgreen};
display_legend = true;
max_val_in_plot = 1;

%%
[x_arrays_iter, quantiles_iter] =  compute_and_plot_LIBSVM_quantiles_in_logscale(num_iter_array, arrsIter, ...
    num_methods, line_colors, display_names, ...
    minmax_colors, quant_colors, display_legend, max_val_in_plot);
ylabel('RSE','FontSize',15)
xlabel('Iter','FontSize',15)
txt=title( ['{\tt ijcnn1}',', $m=$ ',num2str(m),', $n=$ ',num2str(n)]);
set(txt, 'Interpreter', 'latex','FontSize',17);

%% plot the CPU time
maxCPU_AmRABK=max(max(CPU_AREABK));
maxCPU_mRABK=max(max(CPU_AmREABK));
maxCPU_RABK=max(max(CPU_REABK));
maxCPU=max(maxCPU_RABK,max(maxCPU_mRABK,maxCPU_AmRABK));

xlabel_i=maxCPU/opts.Max_iter*[1:opts.Max_iter];
num_iter_array=xlabel_i';

%%
[x_arrays_iter, quantiles_iter] =  compute_and_plot_LIBSVM_quantiles_in_logscale(num_iter_array, arrsIter, ...
    num_methods, line_colors, display_names, ...
    minmax_colors, quant_colors, display_legend, max_val_in_plot);
ylabel('RSE','FontSize',15)
xlabel('CPU time','FontSize',15)
xlim([0,0.5])
txt=title( ['{\tt ijcnn1}',', $m=$ ',num2str(m),', $n=$ ',num2str(n)]);
set(txt, 'Interpreter', 'latex','FontSize',17);


