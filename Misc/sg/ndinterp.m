function Y=ndinterp(out,X)
    if ~iscell(out)
        out={out};
    end
    m=numel(out);
    CYIM=out{1}.CYIM;
    B=out{1}.B;
    B=permute(B,[4,1,2,3]);
    resl=out{1}.resl;
    N=out{1}.N;
    if nargin<2
        nout=prod(resl);
        X=[];
        type=0;
    else
        nout = size(X,1);
        type=1;
    end
    Y=zeros(nout,m);
    for i=1:1:size(CYIM,1)
        c_i=CYIM(i,1);
        %y_i=CYIM(i,2);
        I_i=CYIM(i,3:(2+N));
        M_i=CYIM(i,(3+N:end));     
        %one measurement
        b_i=ones(nout,1);
        for d=1:1:N
            if type==0
                lin = repmat(B(1:resl(d),I_i(d),M_i(d),d)',prod(resl(1:d-1)),1);
                lin = reshape(lin,[],1);
                b_d= repmat(lin,prod(resl)/prod(resl(1:d)),1);
            else
                b_d=B(X(:,d),I_i(d),M_i(d),d);
            end
            b_i=b_i.*b_d;
        end
        for j=1:m
            Y(:,j)=Y(:,j)+c_i*out{j}.CYIM(i,2)*b_i;
        end
    end
end