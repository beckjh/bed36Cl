function histogram_plot(states,resolution_factor,repeats,func,label,confidence)
A = cellfun(func,states);
A = sort(A);
rem = isnan(A);
A(rem) = [];
repeats(rem) = [];
cum_repeats = cumsum(repeats);
sum_repeats = sum(repeats);
cum_weight = cum_repeats/sum_repeats;
upper_quantile = A(min_such_that(cum_weight,@(x) x>(1-(100-confidence)/2/100)));
median = A(min_such_that(cum_weight,@(x) x>=0.5));
lower_quantile = A(min_such_that(cum_weight,@(x) x>=(100-confidence)/2/100));
N_bins = ceil(resolution_factor/4*sum(repeats)^(1/3));
[n,edges,bin] = histcounts(A,N_bins);
W = zeros(size(n));
n_inf = sum(repeats(isinf(A)));
for i = 1:length(repeats)
    if bin(i)>0
        W(bin(i)) = W(bin(i))+repeats(i);
    end
end
if n_inf>0
    warning('%.1f %% of data points for %s were infinite. Plotting only finite data',100*n_inf/sum(repeats),label);
end
W = W/sum(repeats);
h = bar((edges(1:end-1)+edges(2:end))/2,W,'hist');
hold on
h.FaceColor = [150,150,150]/255.0;
ylim_old = ylim;
plot([lower_quantile,lower_quantile],ylim_old,'r');
plot([median,median],ylim_old,'k','linewidth',2);
plot([upper_quantile,upper_quantile],ylim_old,'r');
ylim(ylim_old)
ylabel('Frequency','interpreter','latex')
xlabel(label,'interpreter','latex');
fprintf('%.0f%% quantile, median, and %.0f%% quantile of %s: %.3g, %.3g, %.3g\n',(100-confidence)/2,100-(100-confidence)/2,label,lower_quantile,median,upper_quantile)
end
