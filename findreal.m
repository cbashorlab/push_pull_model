function [pdf,mean_y] = findreal(histogram,xbins,xbin,ybins,ybin,object)
low=1;
high=xbins;

index=binary(xbin,low,high,object);

pdf=histogram(:,index);

mean_y=zeros(ybins,1);

for i=1:ybins
    mean=(ybin(i)+ybin(i+1))/2;
    mean_y(i)=mean;
end

end

