function out = simulateCl36(fault,settings,scenario)
	% Calculate 36Cl concentrations in samples, given a sequence of earthquakes.
	%
	% FAULT is structure with general experiment setup information:
	%   .samples: sample compositions, thicknesses, Cl36 concentrations and standard
	%   deviations
	%   .coll: colluvium composition
	%   .mag_field: magnetic field description
	%   .parameters: alpha, beta, gamma, rho_rock
	%
	% SCENARIO describes uncertain properties of the fault
	%
	%  Chemical production formulae are mostly taken from supplement of
	%   Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
	%   Schimmelpfennig I., Finkel R., Pou K. G.J.Int., 2010
	c = constants();
	par = fault.parameters;
	ppm = ppm_function(fault.samples.chemistry,c.A_k);
	ppmc = ppm_function(fault.coll,c.A_k);
	thick = max(fault.samples.thickness);
	sc_sp = @(x) scaling(x(:,1),x(:,2),thick,x(:,4),x(:,4),par.alpha,par.beta,par.gamma,x(:,3),par.rho_rock);
	sc_mu = @(x) scaling(x(:,1),x(:,2),thick,x(:,4),x(:,4),par.alpha,par.beta,par.gamma,x(:,3),par.rho_rock);
	sc_th = @(x) scaling(x(:,1),x(:,2),thick,repmat(attenuationlengths(ppmc,x(1,3),c,'th'),size(x,1),1),x(:,4),par.alpha,par.beta,par.gamma,x(:,3),par.rho_rock);
	sc_eth = @(x) scaling(x(:,1),x(:,2),thick,repmat(attenuationlengths(ppmc,x(1,3),c,'eth'),size(x,1),1),x(:,4),par.alpha,par.beta,par.gamma,x(:,3),par.rho_rock);
	if nargin == 2
	    out = offline_phase(fault,sc_sp,sc_mu,sc_th,sc_eth,settings);
	    return
    elseif isfield(fault,'offline_data')
        sc_sp = @(x) tensor_evaluation(x,fault.offline_data.sc_sp{1},fault.offline_data.sc_sp{2});
        sc_mu = @(x) tensor_evaluation(x,fault.offline_data.sc_mu{1},fault.offline_data.sc_mu{2});
        sc_th = @(x) tensor_evaluation(x,fault.offline_data.sc_th{1},fault.offline_data.sc_th{2});
        sc_eth = @(x) tensor_evaluation(x,fault.offline_data.sc_eth{1},fault.offline_data.sc_eth{2});
	end
	h = fault.samples.h;
	age = [settings.T_min;scenario.times;0]';
	slip = [0;scenario.jumps;0]';
	N_samples = length(h);
	N = zeros(N_samples,1);
	i_erosion = max_such_that(age,@(x) x<scenario.T_init);
	assert(abs(sum(slip(i_erosion+1:end))-fault.parameters.H_sc)<1);
	depths = sum(slip)-h;%below colluvium along beta plane
	L_th_rock = attenuationlengths(ppm,scenario.rho_rock,c,'th');
	L_eth_rock = attenuationlengths(ppm,scenario.rho_rock,c,'eth');
	[p_cosmo_sp,p_cosmo_mu,p_cosmo_th,p_cosmo_eth,p_rad]=productionCl36(fault.mag_field,scenario,ppm,c);
	heights = sum(slip(1:i_erosion))+par.H_sc-h;%below gamma plane, along beta plane
	for i = 1:(length(age)-1)
	    T = age(i+1)-age(i);
	    depths = depths - slip(i);
	    if i<=i_erosion
		heights = heights - slip(i);
	    end
	    p_sp = p_cosmo_sp.*eval_scaling(sc_sp,par,depths,heights,scenario.rho_coll,scenario.Lambda_sp);
	    p_mu = p_cosmo_mu.*eval_scaling(sc_mu,par,depths,heights,scenario.rho_coll,scenario.Lambda_mu);
	    p_th = p_cosmo_th.*eval_scaling(sc_th,par,depths,heights,scenario.rho_coll,L_th_rock);
	    p_eth = p_cosmo_eth.*eval_scaling(sc_eth,par,depths,heights,scenario.rho_coll,L_eth_rock);
	    p_tot = p_sp+p_mu+p_th+p_eth+p_rad;
	    N = (1./c.lambda36).*exp(-c.lambda36.*T).*(p_tot.*(exp(c.lambda36.*T)-1)+c.lambda36.*N);
	end
	out = N;
end
function y = eval_scaling(sc_function,par,depths,heights,rho_colls,Lambdas)
    N_samples = numel(depths);
    if numel(rho_colls)<N_samples
        rho_colls = repmat(rho_colls,N_samples,1);
    end
    if numel(Lambdas)<N_samples
        Lambdas = repmat(Lambdas,N_samples,1);
    end
    x = [max(0,sind(par.beta-par.alpha)*depths),heights,rho_colls,Lambdas];
    y = sc_function(x);
end
