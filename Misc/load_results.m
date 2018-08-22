function results = load_results(input)
    if ischar(input)
        input_file = ['CaseStudies/',input,'/results.mat'];    
        fprintf('Loading results from %s ...',input_file);
        tmp = load(input_file);
        results = tmp.results;
        fprintf(' done\n');
    else
        if length(fieldnames(input)) == 1
            results = input.results;
        else
            results = input;
        end
    end
end