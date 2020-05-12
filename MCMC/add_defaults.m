function [fault,settings]=add_defaults(case_study)
    %%%% The values in this block can be changed in the case study file%%%%
    prior.recurrence_mean_min = 200;      % Average interevent times are drawn from Inverse Gamma distribution
    prior.recurrence_mean_max = 2000;     % with parameters recurrence_mean and recurrence_alpha;
    prior.recurrence_alpha_min = 1;       % these parameters are inferred within the ranges given here.
    prior.recurrence_alpha_max = 10;      % ... 
    prior.switch_distance_min = 3e3;      % The distance between different long-scale slip-activity regions is drawn from an exponential distribution with mean switch_distance.
    prior.switch_distance_max = 3e4;      % This parameter is inferred within the range given here. 
    prior.tau_min = 0.5;                  % Minimal short-scale variability of inter-event times
    prior.tau_max = 1.5;                  % Maximal ...
    prior.d_min = 10;                     % Minimal displacement (along fault plane)
    prior.d_max = 300;                    % Maximal displacement (along fault plane)
    prior.T_init_min = -20000;            % Lower bound on demise of LGM 
    prior.T_init_max = -12000;            % Upper bound on demise of LGM
    prior.p_zero_freq_max = 0;            % Likelihood of a long-scale activity region with zero activity 
    prior.no_more_slips = 0;              % Time after which no more earthquakes occured
    prior.previous_earthquakes = [];      % 2xN array of known earthquake times and displacement sizes
    settings.T_min = -30000;              % Start of simulations
    settings.debug = false;               % Display debug information during runs
    settings.group_size = 2;              % Number of independent Markov Chains. At least two are required to perform the Gelman-Rubin convergence test 
    settings.pt_levels = gamma(linspace(2,2+19/6,20)); % Parallel tempering levels (each of the independent Markov Chains is run with a number of levels given by the length of this array) 
    settings.dT = 25;                     % Time discretization in simulations
    settings.rho_max = 0.1;               % Model error is inferred between 0 and rho_max
    settings.correlation_length_max = 0;  % Correlation length of model error
    settings.modelscarp =  @simulateCl36; % Function that simulates Cl36 concentrations    
    settings.L = 6;                       % Accuracy of sparse grid interpolation in offline phase    
    settings.thinning = 5;                % One out of how many samples should be saved
    parameters.Psi_sp = @(~)NormalProposal(48.8,1.7); % Production rate of Cl36 by spallation of Ca40  (atoms/g/year), see Gosse & Phillips (2001) for this and the following parameters
    parameters.Psi_mu = @(~)NormalProposal(190,19);   % Production rate of Cl36 by muon capture (atoms/g/year)
    parameters.Lambda_sp = @(~) UniformProposal(180,220);   % Attenuation length for spallation (fast neutron production)
    parameters.Lambda_mu = @(~) UniformProposal(1300,1700); % Attenuation length muon capture
    truth = struct();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case_study_function = str2func(case_study);
    [parameters,prior,settings,truth] = case_study_function(parameters,prior,settings,truth);
    rock = load(['CaseStudies/',case_study,'/rock.txt']);
    colluvium = load(['CaseStudies/',case_study,'/colluvium.txt']);
    magfield = load(['CaseStudies/',case_study,'/magfield.txt']);
    settings.cutoff = 3e4;%How long to project into future to find next earthquake
    settings.T_max = 0;
    try
        ps = parallel.Settings;
        ps.Pool.AutoCreate = ~settings.debug;
    catch e
        msg = 'Error occurred in add_defaults, parallelization not activated (Possible cause: missing the required add-on Parallel Computing Toolbox)';
        warning(msg)
    end
    settings.T_min = min(prior.T_init_min,settings.T_min);
    settings.N_chains = settings.group_size*length(settings.pt_levels);
    prior.lambda_switches_max = (settings.T_max-settings.T_min)/prior.switch_distance_min;
    prior.lambda_switches_min = (settings.T_max-settings.T_min)/prior.switch_distance_max;
    prior.freq_mean_min = 1/prior.recurrence_mean_max;
    prior.freq_mean_max = 1/prior.recurrence_mean_min;
    prior.freq_alpha_min = prior.recurrence_alpha_min;
    prior.freq_alpha_max = prior.recurrence_alpha_max;
    prior.previous_earthquakes = sortrows(prior.previous_earthquakes);
    samples.h = rock(:,end-3)-parameters.H_tr;
    samples.y = rock(:,end-1);
    samples.sigma = rock(:,end);
    samples.thickness = rock(:,end-2);
    samples.chemistry = rock(:,1:62);
    [~,ids] = sort(samples.h);
    samples.h = samples.h(ids);
    samples.chemistry = samples.chemistry(ids,:);
    samples.thickness = samples.thickness(ids,:);
    samples.y = samples.y(ids);
    samples.sigma = samples.sigma(ids);
    fault.coll = colluvium;
    fault.mag_field = magfield;
    parameters.H_sc = round(parameters.H_sc);
    parameters.H_tr = round(parameters.H_tr);
    fault.parameters=parameters;
    fault.prior = prior;
    fault.samples = samples;
    fault.truth = truth;
    if settings.rho_max>0
        fault.parameters.rho = @(~) UniformProposal(0,settings.rho_max+eps);
    else
        fault.parameters.rho = 0;
    end
    if settings.correlation_length_max > 0
        fault.parameters.correlation_length = @(~) UniformProposal(0,settings.correlation_length_max+eps);
    else
        fault.parameters.correlation_length = 0;
    end
end

