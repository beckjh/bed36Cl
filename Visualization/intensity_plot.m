function intensity_plot(states,repeats,x_init,x_end,resolution_factor,plot_time)
    N_bins_x = ceil(resolution_factor*8*sum(repeats)^(1/4)); 
    intensity = zeros(N_bins_x,1);
    x_total = x_end-x_init;
    delta_x = x_total/N_bins_x;
    x = [x_init,x_init+delta_x/2+delta_x*(0:N_bins_x-1),x_end];
    for n = 1:length(states)
        state = states{n};
        times = state.times(2:end);
        jumps = state.jumps(2:end);
        for m = 1:length(times)
            n_bin = floor((times(m)-x_init)/delta_x)+1;
            n_bin = max(1,min(N_bins_x,n_bin));
            intensity(n_bin) = intensity(n_bin)+repeats(n)*jumps(m);
        end
    end
    intensity = intensity/sum(repeats)/(x_total/N_bins_x);
    plott(intensity,[0,0,0],2)
    function plott(l,color,lw)
        y = [l(1);l;l(end)];
        h1 = plot(x/1000,y);
        h1.Color = color;
        h1.LineWidth = lw;
        xlabel('Time (kyr)','interpreter','latex')
        ylabel('Intensity (cm/yr)','interpreter','latex')
        xlim([plot_time(1)/1000.0,plot_time(2)/1000.0]);
        hold on
    end
end