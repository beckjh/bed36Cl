function [fault,settings,output_file,results,final_states] = load_case_study(case_study)
    output_directory = ['CaseStudies/',case_study];
    [fault,settings] = add_defaults(case_study);
    output_file = [output_directory,'/results.mat'];
    if exist(output_file,'file') == 2
        fprintf('Loading previous Markov chains from %s ...',output_file);
        tmp = load(output_file);
        results = tmp.results;
        final_states = tmp.final_states;
        if ~isequal(reduce_fault(results.fault),reduce_fault(fault))||...
                ~isequal(reduce_settings(results.settings),reduce_settings(settings))
            error('Previous results were created with different configuration')
        else
            fault = results.fault;
        end
        fprintf(' done\n')
    else
        if exist(output_directory,'dir') ~= 7
            fprintf('Creating output directory %s\n',output_directory);
            mkdir(output_directory)
        end
        results = [];
        final_states = [];
    end
    offline_file = [output_directory,'/offline_data.mat'];
    if exist(offline_file,'file') == 2
        fprintf('Loading offline data from %s ...',offline_file);
        tmp = load(offline_file);
        fault.offline_data = tmp.offline_data;
        fprintf(' done\n');
    else
        fprintf('Performing offline calculations ...')
        try
            offline_data = settings.modelscarp(fault,settings);
            offline_file = [output_directory,'/offline_data.mat'];
            save(offline_file,'offline_data','-v7.3');
            fault.offline_data = offline_data;
        catch
            fault.offline_data = [];
        end
        fprintf(' done\n')
    end
end
function rf = reduce_fault(fault)
    rf.prior.no_more_slips = fault.prior.no_more_slips;
    rf.parameters.H_sc = fault.parameters.H_sc;
    rf.parameters.H_tr = fault.parameters.H_tr;
    rf.samples = fault.samples;
end
function rs = reduce_settings(settings)
    rs = settings;
    rs.model = [];
    rs.likelihood = [];
    rs.likelihood_from_simulation = [];
end
