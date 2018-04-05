function v = w_quantile(x,w,q)
    x = sort(x,1);
    s = cumsum(w)/sum(w);
    ids = zeros(size(q));
    for j = 1:length(ids)
        ids(j) = min_such_that(s,@(x) x>q(j));
    end
    v = x(ids,:);
end
