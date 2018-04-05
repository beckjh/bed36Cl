function out=sgtensor(out,resl,ex)
    if nargin==2
        ex=zeros(out.N,1);
    end
    X=zeros(max(resl),4);
    for d=1:1:out.N
        if(ex(d))>0
            Yb=exp(-ex(d)*out.ab(1,d));
            Ye=exp(-ex(d)*out.ab(2,d));
            X(1:resl(d),d)=log(1-linspace(Yb,Ye,resl(d)))/(-ex(d));
        end
        X(1:resl(d),d)=linspace(out.ab(1,d),out.ab(2,d),resl(d));
    end
    
    N=out.N;
    w=out.w;
    knots_lvls=out.knots_lvls;
    %w=max(out.I(:));
    B=zeros(w+1,i2m(w+1),N,size(X,1));
    for d=1:N
    for L=1:1:max(out.I(:,d)+1)
        M=i2m(L);
        ma=max(knots_lvls(L,1:M,d));
        mi=min(knots_lvls(L,1:M,d));
        knots_lvls(L,1:M,d)=(knots_lvls(L,1:M,d)-mi)/(ma-mi);
        Xt=(X(1:resl(d),d)-mi)/(ma-mi);
        for j=1:1:M
                B(L,j,d,1:resl(d))=lagrangeb(knots_lvls(L,1:M,d),Xt,knots_lvls(L,j,d));
        end
    end
    end
    out.B=B;
    out.resl=resl;
end