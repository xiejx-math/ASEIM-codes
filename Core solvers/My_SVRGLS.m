function [x,Out]=My_SVRGLS(A,b,opts)

% SVRG for solving least square
%              min 1/2 \|Ax-b\|^2_2
%
%
%Input: the coefficent matrix A, the vector b, size of the block ell,
%and opts
%opts.initial: the initial vector x^0
%opts.TOL: the stopping rule
%.....
%
%Output: the approximate solution x and Out
% Out.error: the relative iterative error \|x^k-x^*\|^2/\|x^k\|^2
% Out.iter: the total number of iteration
% ....
%
% Based on the manuscript: On adaptive stochastic extended iterative methods
% for solving least squares, Yun Zeng, Deren Han, Yansheng Su,
% and Jiaxin Xie, https://arxiv.org/abs/2405.19044
%
% Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn
%

tic
[m,n]=size(A);

%% setting some parameters
flag=exist('opts');
%%%% setting the max iteration
if (flag && isfield(opts,'Max_iter'))
    Max_iter=opts.Max_iter;
else
    Max_iter=20000;
end
%%%% setting the tolerance
if (flag && isfield(opts,'TOL'))
    TOL=opts.TOL;
else
    TOL=10^-12;
end

%%%% setting the tolerance for checking the values of S_k(Ax_k-b)
if (flag && isfield(opts,'TOL1'))
    TOL1=opts.TOL1;
else
    TOL1=10^-20;
end

%%%% setting the initial point
if (flag && isfield(opts,'initial'))
    x=opts.initial;
else
    x=zeros(n,1);
end


%%%% determining what to use as the stopping rule
if (flag && isfield(opts,'xstar'))
    xstar=opts.xstar;
    if m>=n
        normxstar=norm(xstar)^2;
        error1=norm(xstar-x)^2/normxstar;
        strategy=1;
    else
        strategy=0;
    end
else
    strategy=0;
end

if (flag && isfield(opts,'strategy'))
    strategy=opts.strategy;
    normxstar=norm(xstar)^2;
end

if ~strategy
    normb=norm(b)^2+1;
    error1=norm(A*x-b)^2/normb;
end

%%
RSE(1)=error1;

para1=max(sum(A.^2,2));
stepsize_alpha=0.1/para1;


%% executing SAGA
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;
    x0=x;
    u=A'*(A*x0-b)/m;

    for i=1:2*m
        %%% select an index
        l=randperm(m,1);
        g=(A(l,:)*(x-x0))*A(l,:)'+u;
        x=x-stepsize_alpha*g;
    end

    %%% stopping rule
    if strategy
        error1=sum((x-xstar).^2)/normxstar;
        RSE(iter+1)=error1;
        if error1<TOL  || iter>=Max_iter
            stopc=1;
        end
    else
        %%%% Note that we do not us this stopping rule during our test
        error1=norm(A*x-b)^2/normb;
        RSE(iter+1)=error1;
        if  error1<TOL || iter>=Max_iter
            stopc=1;
        end
    end

    times(iter+1)=times(iter)+toc;
end
%% setting Output
Out.error=RSE;
Out.iter=iter;
Out.times=times;
end

