function recent_earthquakes_plot(states,repeats,plot_time,label,resolution_factor)
    A = cellfun(@(x) length(x.times(x.times>plot_time(1))),states);
    n_max = max(A);
    n_max = min(4,n_max);
    cl = colormap('lines');
    total_mass_in_window = zeros(n_max,1);
    handles = {};
    legends = {};
    for j = n_max:-1:1
        A = cellfun(@(x) get_earthquake(x,j),states);
        isn = isnan(A)|(A<plot_time(1)/1000)|(A>plot_time(2)/1000);
        A(isn) = [];
        mod_repeats = repeats;
        mod_repeats(isn) = [];
        N_bins = ceil(resolution_factor*(plot_time(2)-plot_time(1))/(1000*std(A))*sum(mod_repeats)^(1/3)/30);
        N_bins = min(1000,N_bins);
        delta_x = (plot_time(2)-plot_time(1))/N_bins/1000;
        if delta_x > 0 && N_bins > 0     
            edges = plot_time(1)/1000+delta_x*(0:N_bins);
            edges(end) = plot_time(2)/1000;%To fix round off error which puts earthquakes at 0 outside of last edge
            [n,edges,bin] = histcounts(A,edges);
            W=zeros(size(n));
            for i = 1:length(mod_repeats)
                if bin(i) >= 1 && bin(i) <= length(W);
                    W(bin(i)) = W(bin(i))+mod_repeats(i);
                end
            end
            total_mass_in_window(j) = sum(W)/sum(repeats);
            W=W/delta_x/sum(repeats);
            hold on
            y = [W(1),W,W(end)];
            x = [edges(1),edges(2:end)-delta_x/2,edges(end)];
            h = area(x,y);
            handles{end+1} = h;
            legends{end+1} = sprintf('%s most recent event (Total mass in window: %.2f)',num2ordinal_fix(j),total_mass_in_window(j));
            h.FaceColor=cl(j,:);
            try
                h.FaceAlpha = 0.5;
            end
            xlim([plot_time(1)/1000 plot_time(2)/1000]);
        end
    end
    legend([handles{:}],legends);
    xlabel(label,'interpreter','latex');
    box on
end
function f = num2ordinal_fix(x)
    if x==1
        f = 'The';
    else
        f = num2ordinal(x);
    end
end
function a=get_earthquake(state,n)
    if length(state.times) >= n
        a = state.times(end+1-n)/1000.0;
    else
        a = nan;
    end
end