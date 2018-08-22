function slip_history_plot(states,repeats,fault,settings,plot_time,i_max,confidence)
	N_bins = (settings.T_max-settings.T_min)/settings.dT+1;
    times = linspace(settings.T_min,settings.T_max,N_bins)';  
    i_start = min_such_that(times,@(x) x>=plot_time(1));
    N_bins = N_bins+1-i_start;
    times = times(i_start:end);
    leg = {};
    handles = {};
    slip_history_all = zeros(N_bins,0);
    for i = 1:length(states)
        state = strip_pre_erosion(states{i});
        cum_slip = [0;cumsum(state.jumps)];
        ev = zeros(length(times),1);
        for j=1:length(state.times)
            ev = ev+(times>=state.times(j));
        end
        slip_history_all(:,i) = cum_slip(ev+1);
        if i == i_max 
            slip_history_ls = slip_history_all(:,i_max);
            x_values = zeros(1,2*(length(state.times)+1));
            x_values(1:2:end) = [settings.T_min,state.times'];
            x_values(2:2:end) = [state.times',settings.T_max];
            y_values = zeros(1,2*(length(state.times)+1));
            y_values(1:2:end) = [0,cumsum(state.jumps)'];
            y_values(2:2:end) = [0,cumsum(state.jumps)'];
            h = plot(x_values/1000,y_values,'--');
            h_max = h;
            h.Color = [0,0,0];
            h.LineWidth = 2;
        end
        if i == 1
            hold on
        end
    end
    if nargin > 2
        slip_history_quantiles = w_quantile(slip_history_all',repeats,...
            [(confidence+(100-confidence)/2)/100,0.5,(100-confidence)/2/100])';
        h = plot(times/1000.0,slip_history_quantiles(:,2));
        h.Color = [0,0,0];
        h.LineWidth = 2;
        leg{end+1} = 'Medians';
        handles{end+1} = h;
        h = plot(times/1000.0,slip_history_quantiles(:,1));
        h.Color = [1,0,0];
        leg{end+1} = sprintf('$%.0f %s$-confidence bands',confidence,'\%');
        handles{end+1} = h;
        h = plot(times/1000.0,slip_history_quantiles(:,3));
        h.Color = [1,0,0];
    end
    if isfield(fault.truth,'times') && isstruct(fault.truth)
        state = strip_pre_erosion(fault.truth);
        cum_slip = [0;cumsum(state.jumps)];
        slip_history_truth = cum_slip(sum(bsxfun(@ge,times,state.times'),2)+1);
        h = plot(times/1000.0,slip_history_truth,'b-');
        h.LineWidth = 2;
        leg{end+1} = 'Truth';
        handles{end+1} = h;
        slip_history_mean = mean(slip_history_all,2);
        fprintf('L2 distance of mean to truth: %f\n',norm(slip_history_truth-slip_history_mean,2)/power(N_bins,1/2))
        fprintf('L2 distance of LS to truth: %f\n', norm(slip_history_truth-slip_history_ls,2)/power(N_bins,1/2))
        fprintf('L2 distance of median to truth: %f\n', norm(slip_history_truth-slip_history_quantiles(:,2),2)/power(N_bins,1/2))
    end
    handles{end+1} = h_max;
    leg{end+1} = 'Least squares solution';
    xlabel('Time (kyr)','interpreter','latex')
    ylabel('Accumulated displacement (cm)','interpreter','latex')
    legend([handles{:}],leg,'Location','NorthWest','interpreter','latex');
    axis([plot_time(1)/1000.0,plot_time(2)/1000.0,0,fault.parameters.H_sc]);
    box on
end