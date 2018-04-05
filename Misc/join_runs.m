function results = join_runs(results1,results2)
    if ~isempty(results1)
        N_chains = length(results1.state_history);
        results.state_history = cell(N_chains,1);
        results.repeat_history = cell(N_chains,1);
        for j = 1:N_chains
            results.state_history{j} = [results1.state_history{j},results2.state_history{j}];
            results.repeat_history{j} = [results1.repeat_history{j},results2.repeat_history{j}];
            results.n_proposed{j} = results1.n_proposed{j} + results2.n_proposed{j};
            results.n_accepted{j} = results1.n_accepted{j} + results2.n_accepted{j};
        end
        results.state_container = containers.Map([results1.state_container.keys results2.state_container.keys],[results1.state_container.values results2.state_container.values]);
        results.fault = results2.fault;
        results.settings = results2.settings;
    else
        results = results2;
    end
    results.fault = rmfield(results.fault,'offline_data');
end