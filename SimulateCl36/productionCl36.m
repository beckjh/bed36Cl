function [P_cosmo_sp,P_cosmo_mu,P_cosmo_th,P_cosmo_eth,P_rad] = productionCl36(mag_field,scenario,ppm1,c)
% Calculate the production of 36Cl in each sample of a profile per unit
% influx. Most equations here are taken from Schlagenhauf et al., 2010.
Psi_Cl36_Ca_0 = scenario.Psi_sp; % spallation production rate of 36Cl at surface from 40Ca
Psi_mu_0 = scenario.Psi_mu;
EL_f = sum(mag_field(:,3).*mag_field(:,2))/sum(mag_field(:,2));
EL_mu = sum(mag_field(:,4).*mag_field(:,2))/sum(mag_field(:,2));
f_c_Ca = (c.Num_k(49)*ppm1(:,62)*1e-6./c.A_k(49))./(sum(bsxfun(@rdivide,bsxfun(@times,c.Num_k,ppm1(:,1:61)),c.A_k),2)*1e-6) ;% Schlagenhauf et al divide by ppm of colluvium
f_c_K = (c.Num_k(51)*ppm1(:,51)*1e-6./c.A_k(51))./(sum(bsxfun(@rdivide,bsxfun(@times,c.Num_k,ppm1(:,1:61)),c.A_k),2)*1e-6) ; % same as above
X = (sum(bsxfun(@times,ppm1(:,1:61),c.S_i.*c.Y_U_n),2))./(sum(bsxfun(@times,c.S_i,ppm1(:,1:61)),2)) ;
Y = (sum(bsxfun(@times,ppm1(:,1:61),c.S_i.*c.Y_Th_n),2))./(sum(bsxfun(@times,c.S_i,ppm1(:,1:61)),2));
C_Ca = ppm1(:,62)*1e-6 ; % Mass concentration of Ca (g of Ca per g of rock) % from ICP
P_sp_Ca = Psi_Cl36_Ca_0*C_Ca ; % result unscaled 36Cl production by spallation of 40Ca (atoms 36Cl g-1 yr-1)(at of Cl36 /g of K per yr) [162 � 24 Evans et al. 1997]
C_K = ppm1(:,51)*1e-6 ; % Mass concentration of K (g of K per g of rock)
P_sp_K = c.Psi_Cl36_K_0*C_K ; % result unscaled 36Cl production by spallation of 39K (atoms 36Cl g-1 yr-1)
Psi_Cl36_Ti_0 = 13 ; % Spallation production rate at surface of Ti(at of Cl36 /g of Ti per yr) [13 � 3 Fink et al. 2000]
C_Ti = ppm1(:,52)*1e-6 ; % Mass concentration of Ti (g of Ti per g of rock)
P_sp_Ti = Psi_Cl36_Ti_0*C_Ti ; % result unscaled 36Cl production by spallation of Ti (atoms 36Cl g-1 yr-1)
C_Fe = ppm1(:,46)*1e-6 ; % Mass concentration of Fe (g of Fe per g of rock)
P_sp_Fe = c.Psi_Cl36_Fe_0*C_Fe ; % result unscaled 36Cl production by spallation of Fe (atoms 36Cl g-1 yr-1)
P_sp = (P_sp_Ca + P_sp_K + P_sp_Ti + P_sp_Fe); % Unscaled Spallation production rate (atoms 36Cl g-1 yr-1)
Y_Sigma_Ca = f_c_Ca.*c.f_i_Ca.*c.f_d_Ca.*c.f_n_Ca ; % 36Cl production per stopped muon 
Y_Sigma_K = f_c_K.*c.f_i_K.*c.f_d_K.*c.f_n_K ; % 36Cl production per stopped muon 
Y_Sigma = Y_Sigma_Ca + Y_Sigma_K ;
P_mu = Y_Sigma*Psi_mu_0; % Unscaled slow negative muon production rate (atoms 36Cl g-1 yr-1)
U = ppm1(:,37) ; % Concentration en Uranium (ppm1)
Th = ppm1(:,35) ; % Concentration en Thorium (ppm1)
P_n_alphan = X.*U + Y.*Th ; % alpha,n reactions
P_n_sf = 0.429*U ; % spontaneous fission
N_k = bsxfun(@rdivide,ppm1(:,1:61),c.A_k)*c.Avogadro*1e-6 ; %concentrations in atom/g ROCK
N_k(:,56) = N_k(:,56)./scenario.rho_rock ; % divided by bulk-rock density according to CHLOE for H
A = sum(bsxfun(@times,c.A_k,N_k),2)./sum(N_k,2) ;
Sigma_sc = sum(bsxfun(@times,c.sigma_sc_k,N_k),2)*1e-24 ; % (Eq 3.22, Gosse & Phillips
D_eth = 1./(3*Sigma_sc.*(1 - 2./(3*A))) ; % (Eq 3.21, Gosse & Phillips, 2001) Epithermal neutron diffusion coefficient (g.cm-2)
D_th = D_eth ; % D_th = 2.99
Sigma_th = sum(bsxfun(@times,N_k,c.sigma_th_k),2)*1e-24 ; %
L_th = sqrt(D_th./Sigma_th) ;
B = sum(bsxfun(@times,c.Xi_k.*c.sigma_sc_k,N_k),2)*1e-24 ; % Scattering rate parameter 
I_eff = sum(bsxfun(@times,c.I_a_k,N_k),2)*1e-24 ; % (Eq 3.9, Gosse & Philli
Xi = B./Sigma_sc ;  % Eq 3.19 Goss and Phillips, Average log decrement energy loss per neutron collision
Sigma_eth = Xi.*(I_eff + Sigma_sc) ; % (Eq 3.18, Gosse & Phillips, 2001) Effective epithermal loss cross-section (cm2.g-1)
L_eth = 1./sqrt(3*Sigma_sc.*Sigma_eth); % Epithermal neutron diffusion length (g cm-2)
p_E_th = exp(-I_eff./B);%
f_eth = bsxfun(@rdivide,N_k(:,61)*c.I_a_k(61)*(1e-24),I_eff);
f_th = c.sigma_th_k(61)*N_k(:,61)*1e-24./Sigma_th ;
P_th_r = (P_n_alphan + P_n_sf).*p_E_th ; % total radiogenic thermal neutron production
P_eth_r = (P_n_alphan + P_n_sf).*(1 - p_E_th) ; % total radiogenic epithermal neutron production
P_rad = P_th_r.*f_th + P_eth_r.*f_eth ;
Y_s = sum(bsxfun(@times,c.f_d_k.*c.Y_n.*c.Num_k./c.A_k,ppm1(:,1:61)),2)./sum(bsxfun(@times,ppm1(:,1:61),c.Num_k./c.A_k),2) ; %original with ppmc
f_eth = bsxfun(@rdivide,N_k(:,61)*c.I_a_k(61)*(1e-24),I_eff); % (Eq 3.17, Gosse & Phillips, 2001) Fraction of epith neutrons absorbed by Cl35
p_E_th = exp(-I_eff./B) ; % (Eq 3.8, Gosse & Phillips, 2001) Resonance escape probability of a neutron from the epith energy range in subsurface
R_eth = sqrt(A./c.A_a) ; % (Eq 3.24, Gosse & Phillips, 2001) Ratio of epithermal neutron production in subsurface to that in atm
Lambda_eth = 1./Sigma_eth ; % (Eq 3.18,Gosse & Phillips, 2001) Attenuation length for absorbtion and moderation of epith neutrons flux (g.cm-2)
D_eth_a = 1./(3*c.Sigma_sc_a.*(1 - 2./(3*c.A_a))) ; % (Eq 3.21, Gosse & Phillips, 2001) Epithermal neutron diffusion coefficient in atmosphere (g.cm-2)
phi_star_eth = c.P_f_0.*R_eth./(Sigma_eth - (D_eth./(scenario.Lambda_sp.^2))) ; % Epithermal neutron flux at land/atmosphere interface that would be observed in ss if interface was not present (n cm-2 yr-1)
phi_star_eth_a = c.P_f_0.*c.R_eth_a./(c.Sigma_eth_a - (D_eth_a./(scenario.Lambda_sp.^2))) ; % Epithermal neutron flux at land/atmosphere interface that would be observed in atm if interface was not present (n cm-2 yr-1)
P_n_mu_0 = (Y_s*Psi_mu_0 + 5.8e-6*c.phi_mu_f_0); % Fast muon flux at land surface SLHL, Eq.3.49 Gosse & Phillips, 2001 (n cm-2 yr-1)
R_mu = EL_mu.*P_n_mu_0./(EL_f.*c.P_f_0.*R_eth) ; %Ratio of muon production rate to epithermal neutron production rate
Deltaphi_2star_eth_a = phi_star_eth - D_eth_a.*phi_star_eth_a./D_eth ; % Adjusted difference between hypothetical equilibrium epithermal neutron fluxes in atm and ss (n cm-2 yr-1)
L_eth_a = 1./sqrt(3*c.Sigma_sc_a.*c.Sigma_eth_a); % Epithermal neutron diffusion length in atm (g cm-2)
FDeltaphi_star_eth = ((D_eth_a./L_eth_a).*(phi_star_eth_a - phi_star_eth) - ...
    Deltaphi_2star_eth_a.*(D_eth./scenario.Lambda_sp))./...
    ((D_eth_a./L_eth_a) + (D_eth./L_eth)) ; % EQ. 3.28 Gosse & Phillips, 2001 Difference between phi_star_eth,ss and actual epithermal neutron flux at land surface
phi_eth_total = phi_star_eth+...%.*exp(-e./Lambda_e) + ...
    (1 + R_mu.*R_eth).*FDeltaphi_star_eth+...%.*exp(-e./L_eth) + ...
    R_mu.*phi_star_eth;%.*exp(-e./c.Lambda_mu) ; % Epithermal neutron flux (concentration) (n cm-2 yr-1)
P_eth = (f_eth./Lambda_eth).*phi_eth_total.*(1 - p_E_th) ;

f_th = c.sigma_th_k(61)*N_k(:,61)*1e-24./Sigma_th ; % Eq 3.32 de Gosse and Phillips, 2001, fraction of thermal neutrons absorbed by Cl35
Lambda_th = 1./Sigma_th ; % Eq 3.35 Gosse anf Phillips, 2001, Attenuation length for absorbtion of thermal neutrons flux (g.cm-2)
R_th = p_E_th./c.p_E_th_a ; % Ratio of thermal neutron production in ss to that in atm ; Eq 3.34 Gosse and Phillips, 2001
L_th_a = sqrt(c.D_th_a/c.Sigma_th_a) ; % thermal neutron diffusion length in atm (g cm-2)
Deltaphi_star_eth_a = phi_star_eth - phi_star_eth_a ; %difference in equilibrium epithermal neutron fluxes between atm and ss
FDeltaphi_star_eth_a = (D_eth.*Deltaphi_star_eth_a./L_eth - D_eth.*Deltaphi_2star_eth_a./scenario.Lambda_sp)./ ...
    (D_eth_a ./ L_eth_a + D_eth ./ L_eth );
phi_star_th = (c.p_E_th_a.*R_th.*phi_star_eth)./(Lambda_eth.*(Sigma_th - D_th./(scenario.Lambda_sp.^2))) ;% thermal neutron flux at land/atm interface that would be observed in atm if interface not present (n.cm_2.a-1)
R_prime_mu = (c.p_E_th_a./p_E_th).*R_mu ; % ratio of muon production rate to thermal neutron production rate
JDeltaphi_star_eth = (c.p_E_th_a.*R_th.*FDeltaphi_star_eth)./(Lambda_eth.*(Sigma_th - D_th./(L_eth.^2))) ; % Eq. 3.39 Gosse & Phillips, 2001, Portion of difference between phi_star_eth,ss and actual flux due to epithermal flux profile
JDeltaphi_star_eth_a = (c.p_E_th_a.*c.R_th_a.*FDeltaphi_star_eth_a)./((1./c.Sigma_eth_a).*(c.Sigma_th_a - c.D_th_a./(L_eth_a.^2))) ;% Portion of difference between phi_star_eth,a and actual flux due to epithermal flux profile
phi_star_th_a = (c.p_E_th_a.*c.R_th_a.*phi_star_eth_a)./(1./c.Sigma_eth_a.*(c.Sigma_th_a - c.D_th_a./(scenario.Lambda_sp.^2))) ; % thermal neutron flux at land/atmosphere interface that would be observed in atm if interface was not present (n cm-2 yr-1)
Deltaphi_star_th = phi_star_th_a - phi_star_th ; % difference between hypothetical equilibrium thermal neutron fluxes in atmosphere and ss
JDeltaphi_star_th = (c.D_th_a.*(phi_star_th_a./scenario.Lambda_sp - JDeltaphi_star_eth_a./L_eth_a) - ...
    D_th.*(phi_star_th./scenario.Lambda_sp + JDeltaphi_star_eth./L_eth) + ...
    (c.D_th_a./L_th_a).*(Deltaphi_star_th + JDeltaphi_star_eth_a - JDeltaphi_star_eth))./ ...
    ((D_th./L_th) + (c.D_th_a./L_th_a)) ; % portion of difference between phi_star_th,ss and actual flux due to thermal flux profile
phi_th_total = phi_star_th+...%.*exp(-e./Lambda_e) + ...
    (1 + R_prime_mu).*JDeltaphi_star_eth+...%.*exp(-e./L_eth) + ...
    (1 + R_prime_mu.*R_th).*JDeltaphi_star_th+...%.*exp(-e./L_th) + ...
    R_prime_mu.*phi_star_th;%.*exp(-e./c.Lambda_mu) ; % Thermal neutron flux (n.cm_2.a-1)
P_th = (f_th./Lambda_th).*phi_th_total ; % Result unscaled sample specific 36Cl production rate by capture of thermal neutrons (atoms 36Cl g-1 yr-1)
P_cosmo_sp = EL_f.*P_sp;
P_cosmo_mu = EL_mu.*P_mu ;
P_cosmo_th = EL_f.*c.S_L_th*P_th;
P_cosmo_eth = EL_f.*c.S_L_th.*P_eth;
end