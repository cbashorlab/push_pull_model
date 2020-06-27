function index = binary(arr,low,high,x)
max_value=max(arr);
min_value=min(arr);
length=size(arr,1);
if x>=max_value
    index=length-1;
elseif x<=min_value
    index=1;
else
    mid=floor(((high+low))/2);
    if x>=arr(mid) && x<arr(mid+1)
        index=mid;
    elseif x<arr(mid)
        index=binary(arr,low,mid-1,x);
    elseif x>=arr(mid+1)
        index=binary(arr,mid+1,high,x);
    end
end
end

