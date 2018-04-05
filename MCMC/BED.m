function [results,final_states]=BED(fault,settings,runtime,initial_states)
    time_start = tic;
    model = @(x) settings.modelscarp(fault,settings,x);
    likelihood = @(state,sigma,sigma_factor) likelihood_full(model,state,sigma,sigma_factor,fault);
    likelihood_from_simulation = @(state,sigma,sigma_factor) likelihood_from_simulation_full(state,sigma,sigma_factor,fault);
    N_chains = settings.N_chains;
    state_history = cell(N_chains,1);
    repeat_history = cell(N_chains,1);
    state_container = containers.Map('keyType','uint64','valueType','any');
    n_states = zeros(N_chains,1);
    sampler = cell(N_chains,1);
    n_proposed = cell(N_chains,1);
    n_accepted = cell(N_chains,1);
    sigma_factors = zeros(N_chains,1);
    if ~exist('initial_states','var')
        initial_states = [];   
    end
    for n_chain = 1:N_chains
        state_history{n_chain} = cell(1,1);
        level = ceil(n_chain/settings.group_size);
        sigma_factors(n_chain) = settings.pt_levels(level);
        if ~isempty(initial_states) && ~isempty(initial_states{n_chain})
            sampler{n_chain} = Sampler(fault,initial_states{n_chain},settings);
        else
            sampler{n_chain} = Sampler(fault,[],settings);
            [sampler{n_chain}.state.likelihood,sampler{n_chain}.state.simulation] = likelihood(sampler{n_chain}.state,fault.samples.sigma,sigma_factors(n_chain));
            id = randi(2^53-1);
            sampler{n_chain}.state.id = id;
        end
        state_history{n_chain}{1} = sampler{n_chain}.state.id;
        n_states(n_chain) = 1;
        repeat_history{n_chain}(1) = 1;
        tmp = sampler{n_chain}.n_proposals;
        n_proposed{n_chain} = zeros(tmp+1,1);
        n_accepted{n_chain} = zeros(tmp+1,1);
    end
    if isempty(initial_states) && settings.debug
        fprintf('Found valid initial states for all chains\n');
    end
    for n_chain = 1:settings.group_size%save state of level 1
        store_state = sampler{n_chain}.state;
        if ~settings.debug
            store_state.bm = store_state.bm(end);
        end
        state_container(sampler{n_chain}.state.id) = store_state;
    end
    n = 1;
    while toc(time_start) < runtime*60
        n = n+1;
        if settings.debug; fprintf('Individual updates \n'); end
        new_states = zeros(N_chains,1);
        ids = cell2mat(state_container.keys);
        for n_chain = 1:N_chains
            [proposal,MH_ratio,type] = sampler{n_chain}.propose();
            tic_likeli = tic;
            if MH_ratio>0%avoid unnecessary computations
                [likelihood_proposal,simulation_proposal] = likelihood(proposal,fault.samples.sigma,sigma_factors(n_chain));
            else
                likelihood_proposal = 0;
                simulation_proposal = [];
            end
            if settings.debug; fprintf('Chain %d computed likelihood of %dth proposal (type %d) with %d events in %.0f seconds\n',n_chain,n,type,length(proposal.times),toc(tic_likeli)); end;
            n_proposed{n_chain}(type) = n_proposed{n_chain}(type)+1;
            MH_ratio = MH_ratio*likelihood_proposal/sampler{n_chain}.state.likelihood;
            if unifrnd(0,1) < MH_ratio
                if settings.debug; fprintf('Chain %d accepted proposal (type %d)\n',n_chain,type); end;
                n_states(n_chain) = n_states(n_chain)+1;
                sampler{n_chain} = Sampler(fault,proposal,settings);
                sampler{n_chain}.state.simulation = simulation_proposal;
                sampler{n_chain}.state.likelihood = likelihood_proposal;
                id = randi(2^53-1);
                while ismember(id,ids)
                    id = randi(2^53-1);
                end
                sampler{n_chain}.state.id = id;
                state_history{n_chain}{n_states(n_chain)} = id;
                repeat_history{n_chain}(n_states(n_chain)) = 1;
                new_states(n_chain) = 1;
                n_accepted{n_chain}(type) = n_accepted{n_chain}(type)+1;
                if settings.debug>1
                    figure(n_chain)
                    subplot(2,1,1)
                    sampler{n_chain}.plot_tension(proposal);
                    drawnow
                    distFig()
                end
            else
                if settings.debug > 1
                    figure(n_chain)
                    subplot(2,1,2)
                    sampler{n_chain}.plot_tension(proposal);
                    drawnow
                    distFig()
                end
                repeat_history{n_chain}(n_states(n_chain)) = repeat_history{n_chain}(n_states(n_chain))+1;
            end     
        end
        for n_chain = 1:settings.group_size%save states of level 1
            if new_states(n_chain)
                store_state = sampler{n_chain}.state;
                if ~settings.debug
                    store_state.bm = store_state.bm(end);
                end
                state_container(sampler{n_chain}.state.id) = store_state;
            end
        end
        type = sampler{1}.n_proposals+1;
        if settings.debug; fprintf('Swap updates \n'); end;
        for n_chain = N_chains-settings.group_size:-1:1
            n_chain_2 = n_chain+settings.group_size;
            if n_chain_2 <= N_chains
                n_proposed{n_chain}(type) = n_proposed{n_chain}(type)+1;
                n_proposed{n_chain_2}(type) = n_proposed{n_chain_2}(type)+1;
                likelihood_1_new = likelihood_from_simulation(sampler{n_chain_2}.state,fault.samples.sigma,sigma_factors(n_chain));
                likelihood_2_new = likelihood_from_simulation(sampler{n_chain}.state,fault.samples.sigma,sigma_factors(n_chain_2));
                MH_ratio = (likelihood_1_new*likelihood_2_new)...
                    /(sampler{n_chain}.state.likelihood*sampler{n_chain_2}.state.likelihood);
                if unifrnd(0,1)<MH_ratio || (isnan(MH_ratio) && sampler{n_chain_2}.state.likelihood>0)
                    if isnan(MH_ratio) && settings.debug;fprintf('By zero convention ');end;
                    n_states(n_chain) = n_states(n_chain)+1;
                    n_states(n_chain_2) = n_states(n_chain_2)+1;
                    temp_state = sampler{n_chain}.state;
                    sampler{n_chain} = Sampler(fault,sampler{n_chain_2}.state,settings);
                    sampler{n_chain_2} = Sampler(fault,temp_state,settings);
                    sampler{n_chain}.state.likelihood = likelihood_1_new;
                    sampler{n_chain_2}.state.likelihood = likelihood_2_new;
                    id1 = state_history{n_chain_2}{n_states(n_chain_2)-1};
                    state_history{n_chain}{n_states(n_chain)} = id1;
                    id2 = state_history{n_chain}{n_states(n_chain)-1};
                    state_history{n_chain_2}{n_states(n_chain_2)} = id2;
                    if ismember(n_chain,1:settings.group_size)
                        store_state = sampler{n_chain}.state;
                        if ~settings.debug
                            store_state.bm = store_state.bm(end);
                        end
                        state_container(id1) = store_state;
                    end
                    repeat_history{n_chain}(n_states(n_chain)) = 1;
                    repeat_history{n_chain_2}(n_states(n_chain_2)) = 1;
                    if settings.debug; fprintf('Chain %d accepted proposal (type %d)\n',n_chain,type);end;
                    n_accepted{n_chain}(type) = n_accepted{n_chain}(type)+1;
                    n_accepted{n_chain_2}(type) = n_accepted{n_chain_2}(type)+1;
                    if settings.debug>1
                        figure(n_chain)
                        subplot(2,1,1)
                        sampler{n_chain}.plot_tension(sampler{n_chain}.state);
                        drawnow
                        figure(n_chain_2)
                        subplot(2,1,1)
                        sampler{n_chain_2}.plot_tension(sampler{n_chain_2}.state);
                        drawnow
                        distFig()
                    end
                else
                    repeat_history{n_chain}(n_states(n_chain)) = repeat_history{n_chain}(n_states(n_chain))+1;
                    repeat_history{n_chain_2}(n_states(n_chain_2)) = repeat_history{n_chain_2}(n_states(n_chain_2))+1;
                end      
            end
        end
    end
    acceptance = cell(N_chains,1);
    n_proposals = sampler{1}.n_proposals;
    for n_chain = 1:N_chains
        if settings.debug
            fprintf('Accepted proposals:\n')
            fprintf('Type %d: %d/%d=%f\n',[(1:n_proposals+1);n_accepted{n_chain}';n_proposed{n_chain}';(n_accepted{n_chain}./n_proposed{n_chain})']);
            fprintf('Total acceptance ratio: %f\n',sum(n_accepted{n_chain})/sum(n_proposed{n_chain}));
        end
        acceptance{n_chain} = n_accepted{n_chain}./n_proposed{n_chain};
    end
    results.state_history = state_history(1:settings.group_size);
    results.repeat_history = repeat_history(1:settings.group_size);
    results.state_container = state_container;
    results.fault = fault;
    results.settings = settings;
    results.n_accepted = n_accepted;
    results.n_proposed = n_proposed;
    final_states = cellfun(@(x) x.state,sampler,'UniformOutput',false);
end

function [likelihood,simulation] = likelihood_full(model,state,sigma,sigma_factor,fault)
    simulation = model(state);
    state.simulation = simulation;
    likelihood = likelihood_from_simulation_full(state,sigma,sigma_factor,fault);
end

function likelihood = likelihood_from_simulation_full(state,sigma,sigma_factor,fault)
    y = fault.samples.y;
    rho = state.rho;
    simulation = state.simulation;
    kernel = rho.^2*(simulation*simulation');
    if any(isnan(simulation)|isinf(simulation))||any(isnan(rho)|isinf(rho))
        save('fail.mat','state')
        error('Error')
    end
    if state.correlation_length>0
        kernel = kernel.*exp(-abs(bsxfun(@minus,fault.samples.h,fault.samples.h'))/state.correlation_length);
    else
        kernel = diag(diag(kernel)); 
    end
    sum_kernel = diag(sigma.^2)+kernel;
    [Q,D] = eig(sum_kernel);
    simulation = Q'*simulation;
    y = Q'*y;
    D = diag(D);
    likelihood = prod((sigma./sqrt(D)).^(1./sigma_factor.^2).*exp(-(simulation-y).^2./(2.*sigma_factor.^2.*D)));
    if ~isempty(fault.prior.previous_earthquakes)
        if length(state.jumps)<size(fault.prior.previous_earthquakes,1)
           likelihood = 0;
        else
            likelihood = likelihood.*prod(exp(-(fault.prior.previous_earthquakes(:,1)-state.times(end+1-size(fault.prior.previous_earthquakes,1):end)).^2./(2.*sigma_factor.^2.*fault.prior.previous_earthquakes(:,2).^2)));
        end
    end
end
