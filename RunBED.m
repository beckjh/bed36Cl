function RunBED(case_study,save_interval,terminate) 
    % MCMC algorithm for Bayesian earthquake displacement inference using 
    % Cl36 measurements along normal faults.
    %
    % CASE_STUDY Name of directory in subdirectory `CaseStudies`
    % SAVE_INTERVAL (optional, default 60) Time between successive saves (min)
    % TERMINATE (optional) Terminate after TERMINATE saves
    %
    % The actual MCMC algorithm is in the function `BED`, which in turn
    % uses the proposals in the class `Sampler` 
    % (both are located within `MCMC` subdirectory)
    %
    % Published as supplement of:  Beck, J., Wolfers, S., & Roberts, G. P. (2018).
    % Bayesian earthquake dating and seismic hazard assessment using chlorine-36 measurements
    if nargin<2
        save_interval = 60;
    end
    add_BED_path()
    [fault,settings,output_file,results,final_states] = load_case_study(case_study);
    fprintf('Starting MCMC simulation for case study %s \n',case_study);
    saves = 0;
    while true
        fprintf('Next intermediate results should be available around %s\n',...
            char(datetime('now')+minutes(save_interval+1)));
        [new_results,final_states] = BED(fault,settings,save_interval,final_states);
        thinned_new_results = thin(new_results,settings.thinning);
        results = join_runs(results,thinned_new_results);
        save(output_file,'results','final_states');
        fprintf('Results stored in %s\n',output_file)
        if any(cellfun(@(x) isempty(x),final_states(1:settings.group_size)))
            warning('Need to run longer to visualize results')
        end
        saves = saves+1;
        if nargin>2
            if saves==terminate 
                fprintf('Terminating BED ...')
                break
            end
        end
    end
    fprintf(' done\n')
end

