% this file can produce Figures 3 and 4 in the manuscript

close all;
clear;

%m=100;
n=1000;
ell=300;
run_times=10;
rank=n;
kappa=10;

%%%
test_interval=[2000:2000:20000];
sizeP=length(test_interval);

SVRGtimes=zeros(run_times,sizeP);
REABKtimes=zeros(run_times,sizeP);
AREABKtimes=zeros(run_times,sizeP);
AmREABKtimes=zeros(run_times,sizeP);

REABKiter=zeros(run_times,sizeP);
AREABKiter=zeros(run_times,sizeP);
AmREABKiter=zeros(run_times,sizeP);

for ii=1:sizeP
    m=test_interval(ii);
    for jj=1:run_times

        %%%
        [U,~]=qr(randn(m, rank), 0);
        [V,~]=qr(randn(n, rank), 0);
        D = diag(1+(kappa-1).*rand(rank, 1));
        A=U*D*V';

        %%%% generated the right-hand vector b
        x=randn(n,1);
        r1=randn(m,1);
        r=r1-U*(U'*r1);
        b=A*x+r;

        xLS=V*(D^(-1)*(U'*b));

        clear U V D
        %%% parameter setup
        opts.xstar=xLS;
        % SVRG
        [xSVRG,OutSVRG]=My_SVRGLS(A,b,opts);

        
        % REABK
        [xREABK,OutREABK]=My_REABK(A,b,ell,opts);

        %
        [xAREABK,OutAREABK]=My_AREABK(A,b,ell,opts);

        % AmREABK
        [xAmREABK,OutAmREABK]=My_AmREABK(A,b,ell,opts);

        %%%
        SVRGtimes(jj,ii)=OutSVRG.times(end);
        REABKtimes(jj,ii)=OutREABK.times(end);
        AREABKtimes(jj,ii)=OutAREABK.times(end);
        AmREABKtimes(jj,ii)=OutAmREABK.times(end);
        %%%
        REABKiter(jj,ii)=OutREABK.iter;
        AREABKiter(jj,ii)=OutAREABK.iter;
        AmREABKiter(jj,ii)=OutAmREABK.iter;
    end
    fprintf('Done, m=%d\n',m)
end


xlable=test_interval;


%%
y2=SVRGtimes';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=REABKtimes';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%
y4=AREABKtimes';
miny4=min(y4');
maxy4=max(y4');
y4q25=quantile(y4,0.25,2);
y4q75=quantile(y4,0.75,2);

%%
y5=AmREABKtimes';
miny5=min(y5');
maxy5=max(y5');
y5q25=quantile(y5,0.25,2);
y5q75=quantile(y5,0.75,2);

%%
figure
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'black','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'black','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny4 fliplr(maxy4)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y4q25' fliplr(y4q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny5 fliplr(maxy5)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y5q25' fliplr(y5q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)
p2=semilogy( xlable, median(y2'), 'black', 'LineWidth', 1,...
    'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'blue', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
p4=semilogy( xlable, median(y4'), 'red', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
p5=semilogy( xlable, median(y5'), 'green', 'LineWidth', 1,...
    'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('CPU time','FontSize',15)
xlabel('Number of rows $(m)$','Interpreter', 'latex','FontSize',15)
legend([ p2 p3,p4,p5],{'SVRG','REABK','AREABK','AmREABK'},'Interpreter', 'latex','location', 'best','FontSize',14)
txt=title(['{\tt randn}',', $n=$ ',num2str(n),', $r=$ ',num2str(rank),', $\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex','FontSize',17);


%%%%%%%
xlable=test_interval;

%%
y1=REABKiter';
miny1=min(y1');
maxy1=max(y1');
y1q25=quantile(y1,0.25,2);
y1q75=quantile(y1,0.75,2);

%%
y2=AREABKiter';
miny2=min(y2');
maxy2=max(y2');
y2q25=quantile(y2,0.25,2);
y2q75=quantile(y2,0.75,2);

%%
y3=AmREABKiter';
miny3=min(y3');
maxy3=max(y3');
y3q25=quantile(y3,0.25,2);
y3q75=quantile(y3,0.75,2);

%%

figure
h = fill([xlable  fliplr(xlable)], [miny1 fliplr(maxy1)],'blue','EdgeColor', 'none');
set(h,'facealpha', .05)
hold on
h = fill([xlable  fliplr(xlable)], [y1q25' fliplr(y1q75')],'blue','EdgeColor', 'none');
set(h,'facealpha', .1)
p1=semilogy( xlable, median(y1'), 'blue', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 'o', 'DisplayName', 'RABK');
h = fill([xlable  fliplr(xlable)], [miny2 fliplr(maxy2)],'red','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y2q25' fliplr(y2q75')],'red','EdgeColor', 'none');
set(h,'facealpha', .1)
h = fill([xlable  fliplr(xlable)], [miny3 fliplr(maxy3)],'green','EdgeColor', 'none');
set(h,'facealpha', .05)
h = fill([xlable  fliplr(xlable)], [y3q25' fliplr(y3q75')],'green','EdgeColor', 'none');
set(h,'facealpha', .1)

p2=semilogy( xlable, median(y2'), 'red', 'LineWidth', 1,...
            'LineStyle', '-','Marker', 's', 'DisplayName', 'RABK');
p3=semilogy( xlable, median(y3'), 'green', 'LineWidth', 1,...
            'LineStyle', '-','Marker', '^', 'DisplayName', 'RABK');
ylabel('Number of iterations','FontSize',15)
xlabel('Number of rows','Interpreter', 'latex','FontSize',14)
legend([p1 p2 p3],{'REABK','AREABK','AmREABK'},'Interpreter', 'latex','location', 'best','FontSize',14)
txt=title(['{\tt randn}',', $n=$ ',num2str(n),', $r=$ ',num2str(rank),', $\kappa=$ ',num2str(kappa)]);
set(txt, 'Interpreter', 'latex','FontSize',17);


