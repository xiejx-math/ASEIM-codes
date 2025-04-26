% This Matlab file is used to compare AmREABK, pinv, and lsqminnorm
close all;
clear;

n=20000;
m1=20000;
ell=500;
kappa=10;
rank=n;% full column rank to ensure that x is the minimum Euclidean norm least-squares solution

run_times=20; % average times
opts.Max_iter=60000;
mm=40000:10000:100000;
valuem=[30000 35000 mm];

%%
lsq_cputime=zeros(run_times,length(valuem));
pinv_cputime=zeros(run_times,length(valuem));
AmREABK_cputime3=zeros(run_times,length(valuem));

for ii=1:length(valuem)
    m=valuem(ii);

    for jj=1:run_times
        [U,~]=qr(randn(m, rank), 0);
        [V,~]=qr(randn(n, rank), 0);
        D = diag(1+(kappa-1).*rand(rank, 1));
        A=U*D*V';
        %% generated the right-hand vector b
        x=randn(n,1);
        r1=randn(m,1);
        r=r1-U*(U'*r1);
        b=A*x+r;
        clear U V D
        fprintf('Done\n')

        %% lsqminnorm
        tic
        xLS=lsqminnorm(A,b);
        time_lsq=toc;

        fprintf('m=%d,n=%d,rank=%d,block size=%d, kappa=%d\n',m,n,rank,ell,kappa)

        fprintf('lsqminnorm: CPU time %2.4f,  error %2.4e\n',time_lsq,norm(x-xLS))

        %% pinv
        tic
        xLSpinv=pinv(A)*b;
        time_pinv=toc;


        fprintf('pinv: CPU time %2.4f,  error %2.4e\n',time_pinv,norm(x-xLSpinv))

        %% parameter setup
        opts.xstar=x;
        opts.TOL=10*10^(-29);

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
        fprintf('AmREABK3: CPU time %2.4f, error %2.4e\n',OutAmREABK3.times(end),norm(x-xAmREABK3))

        %%
        lsq_cputime(jj,ii)=time_lsq;
        pinv_cputime(jj,ii)=time_pinv;
        AmREABK_cputime3(jj,ii)=OutAmREABK3.times(end);





    end
    fprintf('m=%d\n',m)

end
%% plot the CPU time



xlable=valuem;

figure
semilogy(xlable, lsq_cputime)
hold on
semilogy(xlable, pinv_cputime)
semilogy(xlable, AmREABK_cputime3)
legend('{\tt lsqminnorm }','{\tt pinv}','AmREABK','Interpreter', 'latex','location', 'best','FontSize',14)

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
h = fill([xlable  fliplr(xlable)], [miny5 fliplr(maxy5)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y5q25' fliplr(y5q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'magenta', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
p2=semilogy( xlable, median(y2'), 'black', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p5=semilogy( xlable, median(y5'), 'green', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('CPU time','FontSize',15)
xlabel('Number of rows $(m)$','Interpreter', 'latex','FontSize',15)
xlim([15000 100000])
legend([p1 p2,p5],{'{\tt lsqminnorm }','{\tt pinv}','AmREABK'},'Interpreter', 'latex','location', 'best','FontSize',14)
txt=title(['{\tt randn}',', $n=$ ',num2str(n),', $r=$ ',num2str(rank),', $\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex','FontSize',17);
