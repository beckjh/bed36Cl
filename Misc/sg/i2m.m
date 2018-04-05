function [m] = i2m(i)
m=zeros(size(i));
for n=1:1:length(i)
    if i(n)==1,
        m(n)=1;
    else
        m(n)=2*i(n)-1; %2^(i(n)-1)+1;
    end
end

end
