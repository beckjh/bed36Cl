function gelman_rubin_plot(results,settings,p_burn,plot_time)
    % Gelman rubin plot for different amounts of samples (after fixed
    % burn-in)
    %
    % Needs to do its own burn in, to be able to show convergence (In
    % main burn in, separated repeats of a state are accumulated)
    n_p = 1;
    p = (1:n_p)/n_p;
    group_size = results.settings.group_size;
    level = 1;
    slip_history_all = cell(group_size,1);
    next_earthquakes = cell(group_size,1);
    N_bins = round((plot_time(2)-plot_time(1))/settings.dT/2);
    N_repeats = sum(results.repeat_history{1});
    N_burn_in = round(p_burn*N_repeats);
    rep_his = cell(group_size,1);
    for j = 1:group_size
        weights = results.repeat_history{j+(level-1)*group_size};
        start = find(cumsum(weights)>N_burn_in,1);
        burn = N_burn_in-sum(weights(1:start-1));
        rep_his{j}(start) = weights(start)-burn;
        rep_his{j} = results.repeat_history{j+(level-1)*group_size}(start:end);
        state_his = results.state_history{j+(level-1)*group_size}(start:end);
        [slip_history_all{j},times] = slip_history_no_plot(results.state_container,state_his,N_bins,plot_time);
        next_earthquakes{j} = next_earthquakes_no_plot(results.state_container,state_his);
    end
    PSRF_next_earthquake = PSRF(next_earthquakes,rep_his,p);
    PSRF_array = PSRF(slip_history_all,rep_his,p);
    PSRF_array(isnan(PSRF_array))=0;
    semilogy(times/1e3,PSRF_array);
    
    %legend_entries = cellfun(@(x) sprintf('Using %d/%d of available scenarios',x,n_p),num2cell(1:n_p),'uniformoutput',false);%(after %.0f %s burn in) p_burn*100,'%'
    %legend({legend_entries{1:end-1},sprintf('Using all available scenarios')});
    if nargin>3
        xlim([plot_time(1)/1e3 plot_time(2)/1e3])
    end
    xlabel('Time (kyr)','interpreter','latex')
    ylabel('PSRF','interpreter','latex')
    fprintf('Maximal PSRF of cumulative displacements: %f\n',max(PSRF_array(:,end)));
    fprintf('PSRF of next earthquake: %f\n',PSRF_next_earthquake(end));
    box on
end
function next_earthquakes = next_earthquakes_no_plot(state_container,ids)
    next_earthquakes = zeros(1,length(ids));
    for i = 1:length(ids)
       next_earthquakes(i) = state_container(ids{i}).time_after; 
    end
end
function [slip_history_all,times] = slip_history_no_plot(state_container,ids,N_bins,plot_time)
	times = linspace(plot_time(1),plot_time(2),N_bins)';    
    slip_history_all = zeros(N_bins,length(ids));
    state_history_mat = cell2mat(ids);
    slip_history_unq = zeros(N_bins,length(state_history_mat));
    [state_ids,~,u] = unique(state_history_mat);
    for i = 1:length(state_ids)
        state = state_container(state_ids(i));
        state = strip_pre_erosion(state);
        cum_slip = [0;cumsum(state.jumps)];
        ev = zeros(length(times),1);
        for j = 1:length(state.times)
            ev = ev+(times>=state.times(j));
        end
        slip_history_unq(:,i) = cum_slip(ev+1);
    end
    for i = 1:length(ids)
        slip_history_all(:,i) = slip_history_unq(:,u(i));
    end
end