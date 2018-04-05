function [parameters,prior,settings,truth]=FRAT(parameters,prior,settings,truth)
parameters.H_sc = 1570;    % scarp height along fault (cm)
parameters.H_tr = 130;     % trench depth along fault (cm)
parameters.rho_coll = @(~)UniformProposal(1.2,1.8);% colluvial wedge mean density (g cm^{-3})
parameters.alpha = 25;     % colluvium wedge dip (degrees) 
parameters.beta = 53;      % preserved scarp dip (degrees)
parameters.gamma = 28;     % upper eroded scarp dip (degrees)
parameters.rho_rock = 2.7; % scarp rock mean density (g cm^{-3})
prior.d_max = 300;          % maximal displacement (cm)
end

