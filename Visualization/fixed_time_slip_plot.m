function fixed_time_slip_plot(states,repeats,settings)


A=cellfun(@(x) slip_between(x.times,x.jumps,time_1,time_2),states);
[n,edges,bin]=histcounts(A);
W=zeros(size(n));
B=zeros(0,1);
for i=1:length(repeats)
    W(bin(i))=W(bin(i))+repeats(i);
    B=[B;repmat(A(i),repeats(i),1)];
end
W=W/sum(W);
subplot(2,1,1);
times=(edges(1:end-1)+edges(2:end))/2;
h=bar((edges(1:end-1)+edges(2:end))/2,W,'hist');
h.FaceColor=[0.3,0.3,0.3];
title('Posterior for slip in given time range');
xlabel('(cm)','interpreter','latex');
subplot(2,1,2);
ecdf(B);
title('Empirical CDF');
xlabel('(cm)','interpreter','latex');

    function s=slip_between(times,jumps,time_1,time_2)
        i_1=min_such_that(times,@(x) x>=time_1);
        i_2=max_such_that(times,@(x) x<=time_2);
        s=sum(jumps(i_1:i_2));
    end
end