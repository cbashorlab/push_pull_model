function [histogram,xbin,ybin]=noiselevel(stain,FP,xbins,ybins)
noise_1 = log10([stain, FP]);
x=noise_1(:,1);
y=noise_1(:,2);
xmin=min(x);
ymin=min(y);
xmax=max(x);
ymax=max(y);
xbin=zeros(xbins+1,1);
ybin=zeros(ybins+1,1);
xbin(1)=xmin;
xbin(xbins+1)=xmax;
ybin(1)=ymin;
ybin(ybins+1)=ymax;
xdiff=(xmax-xmin)/xbins;
ydiff=(ymax-ymin)/ybins;
for i=1:xbins
    xbin(i+1)=xbin(i)+xdiff;
end

for i=1:ybins
    ybin(i+1)=ybin(i)+ydiff;
end

histogram=zeros(ybins,xbins);

for i=1:xbins
    noise=noise_1;
    if i<xbins
        noise(noise(:,1)<xbin(i),:)=[];
        noise(noise(:,1)>=xbin(i+1),:)=[];
    else
        noise(noise(:,1)<xbin(i),:)=[]; 
        %noise(noise(:,1)>xbin(i+1),:)=[];
    end 
    gfp_bin=noise(:,2);
    %disp(size(gfp_bin));
    [count,bin_edge]=histcounts(gfp_bin, ybin);
    if sum(count)==0
        average=0.5*(xbin(i)+xbin(i+1));
        index=binary(ybin,1,ybins,average);
        count(index)=1;
        frequency=count/sum(count);
        histogram(:,i)=frequency;
    else 
        frequency=count/sum(count);
        histogram(:,i)=frequency;
    end
end
xbin=10.^(xbin);
ybin=10.^(ybin);
end
