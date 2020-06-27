function driver
tic
parpool('local',parcluster('local').NumWorkers);
disp(parcluster('local').NumWorkers);
%% Read and process noise data
matrix_noise = readmatrix('Staining_noise.csv');
myc=matrix_noise(:,7);
GFP=matrix_noise(:,16);

indice=find(myc<=0);
myc(indice)=[];
GFP(indice)=[];
indice=find(GFP<=0);
myc(indice)=[];
GFP(indice)=[];
% scatter(myc,GFP);
% set(gca,'Xscale','log');
% set(gca,'Yscale','log');

data = readmatrix('170RR.csv');
K = data(:,7);
S = data(:,10);
P = data(:,16);
indice=find(K<=0);
K(indice)=[];
S(indice)=[];
P(indice)=[];
indice=find(S<=0);
K(indice)=[];
S(indice)=[];
P(indice)=[];
indice=find(P<=0);
K(indice)=[];
S(indice)=[];
P(indice)=[];

K=K(1:1000);
S=S(1:1000);
P=P(1:1000);

xbins=100;
ybins=100;
[histogram,xbin,ybin]=noiselevel(myc,GFP,xbins,ybins);
%assignin('base','hist',histogram);

va=1:0.1:10;
va=10.^(va);

%combine=cell(size(va));
combine=[];
parfor i=1:size(va,2)
    [mean_alpha,prob_alpha]=multicell(histogram,histogram,histogram,xbins,xbin,ybins,ybin,xbins,xbin,ybins,ybin,xbins,xbin,ybins,ybin,K,S,P,va(i));
    va_store=ones(size(mean_alpha))*va(i);
    %combine{i}=[prob_alpha mean_alpha va_store];
    combine=[combine; [mean_alpha va_store prob_alpha]];
end
%assignin('base','prob1',combine);
x=combine(:,1);
y=combine(:,2);
z=combine(:,3);
fig1=scatter3(x,y,z);
set(gca,'XScale','log');
set(gca,'Yscale','log');
saveas(fig1,'scatterplot.fig');
fig2=figure;
xrange=linspace(log10(min(x)),log10(max(x)),100);
xrange=10.^xrange;
yrange=linspace(log10(min(y)),log10(max(y)),100);
yrange=10.^yrange;
[X,Y]=meshgrid(xrange,yrange);
Z=griddata(x,y,z,X,Y);
contour(X,Y,Z,5);
set(gca,'XScale','log');
set(gca,'Yscale','log');
saveas(fig2,'contour.jpg');
toc
end

function [mean_alpha,prob_alpha]=multicell(histogram_K,histogram_S,histogram_P,kxbins,kxbin,kybins,kybin,sxbins,sxbin,sybins,sybin,pxbins,pxbin,pybins,pybin,K,S,P,va)
%% Alpha Range Definition
alpha_bins=100;
minalpha=0;
maxalpha=20;
alpha_bin=[minalpha:(maxalpha-minalpha)/alpha_bins: maxalpha];
alpha_bin=alpha_bin';
alpha_bin=10.^(alpha_bin);
mean_alpha=zeros(size(alpha_bin,1)-1,1);
for i=1:size(mean_alpha,1)
    mean_alpha(i)=(alpha_bin(i)+alpha_bin(i+1))/2;
end
prob_alpha=zeros(size(mean_alpha));

%% Multi Cell Combination
for i=1:size(K,1)
    total=singlecell(histogram_K,histogram_S,histogram_P,kxbins,kxbin,kybins,kybin,sxbins,sxbin,sybins,sybin,pxbins,pxbin,pybins,pybin,K(i),S(i),P(i),va);
    total(total(:,1)<=0,:)=[];
    total(total(:,2)<=0,:)=[];
    total(:,2)=log10(total(:,2));
    for j=1:size(total,1)
        index=binary(alpha_bin,1,alpha_bins,total(j,1));
        prob_alpha(index)=prob_alpha(index)+total(j,2);
    end
end

prob_alpha=prob_alpha/sum(prob_alpha);
%plot(mean_alpha,prob_alpha);
%assignin('base','prob1',prob_alpha);
%assignin('base','alpha1',mean_alpha);
end

function total=singlecell(histogram_K,histogram_S,histogram_P,kxbins,kxbin,kybins,kybin,sxbins,sxbin,sybins,sybin,pxbins,pxbin,pybins,pybin,k,s,p,va)
[k_pdf,mean_k]=findreal(histogram_K,kxbins,kxbin,kybins,kybin,k);
[s_pdf,mean_s]=findreal(histogram_S,sxbins,sxbin,sybins,sybin,s);
[p_pdf,mean_p]=findreal(histogram_P,pxbins,pxbin,pybins,pybin,p);

[K_1,S_1,P_1]=meshgrid(mean_k,mean_s,mean_p);
c=cat(4,K_1,S_1,P_1);
real_value=reshape(c,[],3);

[k_1,s_1,p_1]=meshgrid(k_pdf,s_pdf,p_pdf);
d=cat(4,k_1,s_1,p_1);
prob_value=reshape(d,[],3);

prob_store=prob_value(:,1).*prob_value(:,2).*prob_value(:,3);
k_store=real_value(:,1);
s_store=real_value(:,2);
p_store=real_value(:,3);
p_store=p_store/va;
alpha_store=((p_store-k_store).*(p_store-s_store))./p_store;
total=[alpha_store,prob_store];
end


