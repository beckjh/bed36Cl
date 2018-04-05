function PSRF=PSRF(states,weights,p)
PSRF = zeros(size(states{1},1),length(p));
for k = 1:length(p)
    M = length(states);
    means = zeros(size(states{1},1),M);
    vars = zeros(size(states{1},1),M);
    N_keep = ceil(p(k)*min(cellfun(@(x) sum(x),weights)));
    for n_chain = 1:M
        weighting = weights{n_chain};
        start = 1;
        keep = max_such_that(cumsum(weighting),@(x) x<N_keep)+1;
        keep = min(length(weighting),keep);
        weighting(keep) = min(N_keep-sum(weighting(1:keep-1)),weighting(keep));
        means(:,n_chain) = sum(bsxfun(@times,weighting(start:keep),states{n_chain}(:,start:keep)),2)/sum(weighting(start:keep));
        if sum(weighting(start:keep))==0
            error('Too few samples')
        end
        vars(:,n_chain) = var(states{n_chain}(:,start:keep),weighting(start:keep)/sum(weighting(start:keep)),2);
    end
    N = N_keep;
    B = N*var(means,0,2);
    W = mean(vars,2);
    V = (N-1)/N*W+(M+1)/(M*N)*B;
    PSRF(:,k) = sqrt(V./W);
end
end