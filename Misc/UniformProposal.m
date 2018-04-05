function [ proposal,MH_ratio ] = UniformProposal(a,b,x,N)
    if nargin<4
        N=1;
    end
    proposal=zeros(N,1);
    for j=1:N
        if nargin<3||isempty(x)||unifrnd(0,1)<0.5
            proposal(j)=unifrnd(a(j),b(j));
        else
            found=false;
            while ~found
                temp=normrnd(x(j),(b(j)-a(j))/4);
                if temp>=a(j) && temp<=b(j)
                    found=true;
                    proposal(j)=temp;
                end
            end
        end
    end
    MH_ratio=1;
end
