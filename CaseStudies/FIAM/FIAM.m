function [parameters,prior,settings,truth]=FIAM(parameters,prior,settings,truth)
parameters.H_sc = 2705;    % scarp height along fault (cm)
parameters.H_tr = 115;     % trench depth along fault (cm)
parameters.rho_coll = @(~)UniformProposal(1.2,1.8);% colluvial wedge mean density (g cm^{-3})
parameters.alpha = 23;     % colluvium wedge dip (degrees) 
parameters.beta = 42;      % preserved scarp dip (degrees)
parameters.gamma = 33;     % upper eroded scarp dip (degrees)
parameters.rho_rock = 2.7; % scarp rock mean density (g cm^{-3})
prior.d_max = 300;         % maximal displacement (cm)
end
