function [x,Out]=My_REABK(A,b,ell,opts)

% randomized extended average block Kaczmarz method with adaptive step-sizes
% for solving linear systems
%              Ax=b
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
% Based on the manuscript:
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
    Max_iter=2000000;
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

if (flag && isfield(opts,'initial'))
    z=opts.initialz;
else
    z=b;
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


%% setting the probability
if (flag && isfield(opts,'probset'))
    probset=opts.probset;
else
    probset=0;
end

if probset
    %% Row
    Aarrs=opts.Aarrs;
    barrs=opts.barrs;
    Row_cumsumpro=opts.Rowcumsumpro;
    %% Colum
    Atarrs=opts.Atarrs;
    Column_cumsumpro=opts.Columncumsumpro;
    %% step-size
    alpha=opts.alpha;
else
    %% Row probability
    normAfro=norm(A,'fro')^2;
    tauR=ceil(m/ell);
    blockAnormfro=zeros(tauR,1);
    %prob=zeros(tau,1);
    beta1_max=0;
    for i=1:tauR
        if i==tauR
            psR=((i-1)*ell+1):1:m;
        else
            psR=((i-1)*ell+1):1:(i*ell);
        end
        ApsR=A(psR,:);
        blockAnormfro(i)=norm(ApsR,'fro')^2;
        Aarrs{i}=ApsR;
        barrs{i}=b(psR);
        beta1_max=max(beta1_max,norm(full(ApsR))^2/blockAnormfro(i));
    end
    probR=blockAnormfro/normAfro;
    Row_cumsumpro=cumsum(probR);
    %% Column probability
    tauC=ceil(n/ell);
    blockAtnormfro=zeros(tauC,1);
    %prob=zeros(tau,1);
    beta2_max=0;
    for i=1:tauC
        if i==tauC
            psC=((i-1)*ell+1):1:n;
        else
            psC=((i-1)*ell+1):1:(i*ell);
        end
        AtpsC=A(:,psC);
        blockAtnormfro(i)=norm(AtpsC,'fro')^2;
        Atarrs{i}=AtpsC;
    end
    probC=blockAtnormfro/normAfro;
    Column_cumsumpro=cumsum(probC);
    beta2_max=max(beta2_max,norm(full(AtpsC))^2/blockAtnormfro(i));
end
beta_max=max(beta1_max,beta2_max);
%alpha=1.75/beta_max;
alpha=1/beta_max;
%% executing the AREBKU method
stopc=0;
iter=0;
times(1)=toc;
while ~stopc
    tic
    iter=iter+1;

    %% update z
    lc=sum(Column_cumsumpro<rand)+1;
    AtindexT=Atarrs{lc};
    TtAtz=AtindexT'*z;

    %normTtAtz=norm(TtAtz)^2;

   % if normTtAtz>TOL1
        dz=AtindexT*TtAtz;
        %norm_dz=norm(dz)^2;
        %muz=normTtAtz/norm_dz;
        z=z-(alpha/blockAtnormfro(lc))*dz;
   % else
    %    z=z;
   % end

    %% update x
    lr=sum(Row_cumsumpro<rand)+1;
    AindexR=Aarrs{lr};
    if lr==tauR
        psR=((lr-1)*ell+1):1:m;
    else
        psR=((lr-1)*ell+1):1:(lr*ell);
    end
    bzindexR=barrs{lr}-z(psR);
    Axbz=AindexR*x-bzindexR;
    %normAxb=norm(Axbz)^2;
   % if normAxb>TOL1
        dk=AindexR'*Axbz;
        %norm_dk=norm(dk)^2;
        %alpha=normAxb/norm_dk;
        %% update x
        x=x-alpha*dk/blockAnormfro(lr);
    %else
    %    x=x;
    %end

    %% stopping rule
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
Out.beta_max=beta_max;
end

