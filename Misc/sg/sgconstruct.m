function out=sgconstruct(N,w,L,func,ab)
    % generate Smolyak grid with level w
    smolyakset=@(i) w*([(i-1),zeros(1,N-numel(i))])';%[5 2 3 2]
    I_delta=multiidx_gen(N,smolyakset,L,1);

    % compute normalized CC knots
    w=max(I_delta(:));
    knots_lvls=zeros(w+1,i2m(w+1),N);
    for l=1:1:(w+1)
        knots=knots_CC(i2m(l),0,1);
        knots=repmat(knots,N,1)';
        knots=repmat(ab(1,:),size(knots,1),1)+knots.*(repmat(ab(2,:),size(knots,1),1)-repmat(ab(1,:),size(knots,1),1));
        knots_lvls(l,1:i2m(l),:)=knots;
    end

    % index set for interpolation with coefficients C
    [I,C]=index_set_I(I_delta);

    % pre-compute function values
    mapObj = containers.Map;
    CYIM=term_values(mapObj,func,I,C,knots_lvls);

    out.CYIM=CYIM;
    out.mapObj=mapObj;
    out.I=I;
    out.C=C;
    out.knots_lvls=knots_lvls;
    out.w=w;
    out.N=N;
    out.ab=ab;
end