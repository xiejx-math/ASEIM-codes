% This file can reproduce Figure 2 in the manuscript
close all;
clear;

ss=7;
m=2^(ss+3);
n=2^ss;
rank=n;
kappa=10;

test_interval=0:1:log2(m);
sizeP=length(test_interval);

run_times=50;% average times

%%
REABKtimes=zeros(run_times,sizeP);
REABKiter=zeros(run_times,sizeP);
AREABKtimes=zeros(run_times,sizeP);
AREABKiter=zeros(run_times,sizeP);
AmREABKtimes=zeros(run_times,sizeP);
AmREABKiter=zeros(run_times,sizeP);

%%
Rho_REABK=zeros(run_times,sizeP);
Rho_AREABK=zeros(run_times,sizeP);
Rho_AmREABK=zeros(run_times,sizeP);
Rho_upperbound=zeros(run_times,sizeP);



for ii=1:sizeP
    ell=2^(test_interval(ii));
    for jj=1:run_times

        [U,~]=qr(randn(m, rank), 0);
        [V,~]=qr(randn(n, rank), 0);
        D = diag(1+(kappa-1).*rand(rank, 1));
        A=U*D*V';
        clear U D V

        %% generated the right-hand vector b
        %A=full(A);
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

         if ell==1
            TOL1=10^(-24);
         end
        %opts.Max_iter=1000000;

        %% AREABK

        [xAREABK,OutAREABK]=My_AREABK(A,b,ell,opts);

        rho_AREABK=compute_rhok(OutAREABK.error);

        %% AmREABK
        [xAmREABK,OutAmREABK]=My_AmREABK(A,b,ell,opts);

         rho_AmREABK=compute_rhok(OutAmREABK.error);

        %% REABK
        [xREABK,OutREABK]=My_REABK(A,b,ell,opts);

        rho_REABK=compute_rhok(OutREABK.error);

        if rank==n
            i_kappa_F=svds(A,1,'smallestnz')^2/norm(A,'fro')^2;
        else
            S=svd(A);
            i_kappa_F=S(rank)^2/norm(A,'fro')^2;
        end

        rho_upperbound=1-i_kappa_F/OutREABK.beta_max;

        %% update and store the results
        %%%%%%
        Rho_REABK(jj,ii)=rho_REABK;
        Rho_AREABK(jj,ii)=rho_AREABK;
        Rho_AmREABK(jj,ii)=rho_AmREABK;
        Rho_upperbound(jj,ii)=rho_upperbound;

        %%
        REABKtimes(jj,ii)=OutREABK.times(end);
        REABKiter(jj,ii)=OutREABK.iter;

        AREABKtimes(jj,ii)=OutAREABK.times(end);
        AREABKiter(jj,ii)=OutAREABK.iter;

        AmREABKtimes(jj,ii)=OutAmREABK.times(end);
        AmREABKiter(jj,ii)=OutAmREABK.iter;
        if OutAmREABK.iter>=2000000
            OutAmREABK.iter
            OutAmREABK.times(end)
            fprintf('Max_iter\n');
            return;
        end

        %%%%%%

    end
    fprintf('ell=%d\n',ell)
end


%% plot the CPU time
xlable=test_interval;

%%
y1=REABKtimes';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=AREABKtimes';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=AmREABKtimes';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('CPU')
xlabel('$\log_2(p)$','Interpreter', 'latex')
%legend('RABK','location', 'best')
legend([p1 p2 p3],{'REABK','AREABK','AmREABK'},'Interpreter', 'latex','location', 'best')
txt=title(['{\tt randn}',',$m=$ ',num2str(m),',$n=$ ',num2str(n),',$r=$ ',num2str(rank),',$\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex');



%% plot the convergence factors

%%
y1=Rho_REABK';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=Rho_AREABK';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=Rho_AmREABK';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%
y4=Rho_upperbound';
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2);
y4q75=quantile(y4,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
p4=semilogy( xlable, median(y4'), 'black', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
%hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
%%
ylabel('Convergence factors')
xlabel('$\log_2(p)$','Interpreter', 'latex')
%legend('RABK','location', 'best')
legend([p4 p1 p2 p3],{'Upper bound','REABK','AREABK','AmREABK'},'Interpreter', 'latex','location', 'southeast')
txt=title(['{\tt randn}',',$m=$ ',num2str(m),',$n=$ ',num2str(n),',$r=$ ',num2str(rank),',$\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex');
%%
axesposition=[0.21,0.17,0.3,0.32];
axes('position',axesposition);
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
p4=semilogy( xlable, median(y4'), 'black', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
%hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
xlim([1,4])


%% plot the full number of iterations

P_matrix=zeros(run_times,sizeP);
for kk=1:run_times
    P_matrix(kk,:)=2.^(test_interval');
end

REABK_value=REABKiter.*P_matrix/m;
AREABK_value=AREABKiter.*P_matrix/m;
AmREABK_value=AmREABKiter.*P_matrix/m;

%%
y1=REABK_value';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=AREABK_value';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=AmREABK_value';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('$k \cdot \frac{p}{m}$','Interpreter', 'latex')
xlabel('$\log_2(p)$','Interpreter', 'latex')
%legend('RABK','location', 'best')
legend([p1 p2 p3],{'REABK','AREABK','AmREABK'},'Interpreter', 'latex','location', 'northeast')
txt=title(['{\tt randn}',',$m=$ ',num2str(m),',$n=$ ',num2str(n),',$r=$ ',num2str(rank),',$\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex');
%%
axesposition=[0.21,0.57,0.3,0.32];
axes('position',axesposition);
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
xlim([1,4])







