% This Matlab file is used to compare AmREABK, pinv, and lsqminnorm
close all;
clear;

n=1000;

m1=1000;


kappa=4;

rank=n;% full column rank to ensure that x is the minimum Euclidean norm least-squares solution.

run_times=50; % average times
opts.Max_iter=40000;
valuem=[m1:1000:10000];

%%
lsq_cputime=zeros(run_times,length(valuem));
pinv_cputime=zeros(run_times,length(valuem));
AmREABK_cputime1=zeros(run_times,length(valuem));
AmREABK_cputime2=zeros(run_times,length(valuem));
AmREABK_cputime3=zeros(run_times,length(valuem));

for ii=1:length(valuem)
    m=valuem(ii);

    for jj=1:run_times


        [U,~]=qr(randn(m, rank), 0);
        [V,~]=qr(randn(n, rank), 0);
        D = diag(1+(kappa-1).*rand(rank, 1));
        A=U*D*V';
        clear U V D

        %% generated the right-hand vector b
        x=randn(n,1);
        Z=null(full(A)');
        r=Z*randn(size(Z,2),1);
        b=A*x+r;

        %% lsqminnorm
        tic
        xLS=lsqminnorm(A,b);
        time_lsq=toc;

        %% pinv
        tic
        xLSpinv=pinv(A)*b;
        time_pinv=toc;

        %% parameter setup
        opts.xstar=x;
        opts.TOL=3*10^(-29);
        %%
        ell=100;
        stop1=0;
        while ~stop1
            [xAmREABK1,OutAmREABK1]=My_AmREABK(A,b,ell,opts);

            if OutAmREABK1.iter>=opts.Max_iter
                OutAmREABK1.iter
                OutAmREABK1.times(end)
                fprintf('Max_iter 1\n');
                %return;
            else
                stop1=1;
            end

        end
        %%
        ell=200;
        stop3=0;
        while ~stop3
            [xAmREABK3,OutAmREABK3]=My_AmREABK(A,b,ell,opts);
            if OutAmREABK3.iter>=opts.Max_iter
                OutAmREABK3.iter
                OutAmREABK3.times(end)
                fprintf('Max_iter 3\n');
            else
                stop3=1;
            end
        end

        %%
        lsq_cputime(jj,ii)=time_lsq;
        pinv_cputime(jj,ii)=time_pinv;
        AmREABK_cputime1(jj,ii)=OutAmREABK1.times(end);
        AmREABK_cputime3(jj,ii)=OutAmREABK3.times(end);

        fprintf('m=%d,n=%d,rank=%d,block size=%d, kappa=%d\n',m,n,rank,ell,kappa)

        fprintf('AmREABK1: CPU time %2.4f, error %2.4e\n',OutAmREABK1.times(end),norm(x-xAmREABK1))
        fprintf('AmREABK3: CPU time %2.4f, error %2.4e\n',OutAmREABK3.times(end),norm(x-xAmREABK3))
        fprintf('lsqminnorm: CPU time %2.4f,  error %2.4e\n',time_lsq,norm(x-xLS))
        fprintf('pinv: CPU time %2.4f,  error %2.4e\n',time_pinv,norm(x-xLSpinv))

    end
    fprintf('m=%d\n',m)

end
%% plot the CPU time

xlable=valuem;

%%
y1=lsq_cputime';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=pinv_cputime';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=AmREABK_cputime1';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);


%%
y5=AmREABK_cputime3';
miny5=min(y5');
maxy5=max(y5');
y5q25=quantile(y5,0.25,2);
y5q75=quantile(y5,0.75,2);

%%
figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'magenta','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'magenta','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny5 fliplr(maxy5)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y5q25' fliplr(y5q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
p2=semilogy( xlable, median(y2'), 'black', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'red', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
p5=semilogy( xlable, median(y5'), 'green', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('CPU')
xlabel('Number of rows $(m)$','Interpreter', 'latex')
legend([p1 p2 p3,p5],{'{\tt lsqminnorm }','{\tt pinv}','AmREABK $p=100$','AmREABK $p=200$'},'Interpreter', 'latex','location', 'best')
txt=title(['$n=$ ',num2str(n),',$\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex');
