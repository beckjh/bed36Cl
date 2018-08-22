function [scenarios,weights] = LoadBED(input,burn,thinning)
    % INPUT Name of case study or structure `results` from results.mat file
    % BURN Amount of burn in (between 0 and 1). Default 0.1
    % THINNING Thinning factor (>=1). Default 1
    %
    % Returns the results of a BED case study in form of a cell array
    % of scenarios (MCMC samples) and an array of weights corresponding 
    % to each scenario. 
    add_BED_path()
    results = load_results(input);
    if nargin>2
        results = thin(results,thinning);
    end
    if nargin<2
        burn = 0.1;
    end
    [scenarios,repeats] = burn_in(results,burn); 
    weights = repeats/sum(repeats);
end
