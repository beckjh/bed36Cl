function [parameters,prior,settings,truth]=Synthetic(parameters,prior,settings,truth2)
parameters.H_sc = 2705;    % scarp height along fault (cm)
parameters.H_tr = 115;     % trench depth along fault (cm)
parameters.rho_coll = @(~)UniformProposal(1.2,1.8);% colluvial wedge mean density (g cm^{-3})
parameters.alpha = 23;     % colluvium wedge dip (degrees) 
parameters.beta = 42;      % preserved scarp dip (degrees)
parameters.gamma = 33;     % upper eroded scarp dip (degrees)
parameters.rho_rock = 2.7; % scarp rock mean density (g cm^{-3})
prior.d_max = 300;         % maximal displacement (cm)

truth.T_init=-19000;
truth.switches = [-30000,-29300,-26430,-24730,-18220,-17000,-13350,-11140,-10880, -6550, -6060,  -720];
truth.times=[-29540,-29060,-28530,-28370,-28030,-27540,-27350,-27270,-27130,-27090,-26700,-26610,-26170,-25860,-25780,-25510,-25140,-24560,-24490,-23740,-23460,-22950,-22210,-22060,-19510,-19420,-19180,-18700,-16910,-16750,-16570,-16300,-16110,-16020,-15870,-15710,-15660,-15560,-15350,-15080,-15030,-14990,-14780,-14660,-14600,-14570,-14500,-14410,-14320,-13760,-13630,-13550,-13470,-12990,-12690,-12360,-11820,-11400,-11110,-10230, -9990, -9810, -9540, -9430, -9260, -9210, -8970, -8520, -8450, -8350, -8150, -8080, -7920, -7680, -7620, -7310, -7180, -7130, -7010, -6920, -6260, -5980, -4480, -3960, -3280, -2910, -1200,  -570,  -320,  -190,   -10]';
truth.jumps=[ 62, 94, 19, 58, 68, 49, 59, 98,100,102, 16, 46, 48, 60, 80, 46, 90,108, 78, 72, 83, 64, 45, 40, 13, 18, 27, 44, 38, 36, 84, 24, 32, 91, 11, 44, 40, 49, 94, 36, 45, 70, 43, 30, 47, 19, 23, 16, 46, 15, 86, 35, 24, 21, 20,100, 17, 22, 14, 65, 20, 91, 14, 67, 23, 39, 24, 83, 42, 20, 12, 12, 67, 24, 36, 93, 22, 20, 55, 71, 65, 99, 10, 22, 53,85, 41, 13, 14, 75, 12]';
truth.tau=1;
end
