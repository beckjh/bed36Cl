function L=attenuationlengths(ppm,rho,c,which)
% Compute L_th and L_eth, for colluvium or rock, depending on choices of
% inputs ppm and rho.
% Most equations here are taken from Schlagenhauf et al. (2010)
N_k = bsxfun(@rdivide,ppm(:,1:61),c.A_k)*c.Avogadro*1e-6 ; %concentrations in atom/g
N_k(:,56) = N_k(:,56)./rho ; % divided by bulk density according to CHLOE for H
A = sum(bsxfun(@times,c.A_k,N_k),2)./sum(N_k,2) ;
Sigma_sc = sum(bsxfun(@times,c.sigma_sc_k,N_k),2)*1e-24 ; % (Eq 3.22, Gosse & Phillips
if strcmp(which,'th')
    D_eth = 1./(3*Sigma_sc.*(1 - 2./(3*A))) ; % (Eq 3.21, Gosse & Phillips, 2001) Epithermal neutron diffusion coefficient (g.cm-2)
    D_th = D_eth ; % D_th = 2.99
    Sigma_th = sum(bsxfun(@times,N_k,c.sigma_th_k),2)*1e-24 ; %
    L_th = sqrt(D_th./Sigma_th) ;
    L = L_th;
else
    B = sum(bsxfun(@times,c.Xi_k.*c.sigma_sc_k,N_k),2)*1e-24 ; % Scattering rate parameter 
    I_eff = sum(bsxfun(@times,c.I_a_k,N_k),2)*1e-24 ; % (Eq 3.9, Gosse & Philli
    Xi = B./Sigma_sc ;  % Eq 3.19 Goss and Phillips, Average log decrement energy loss per neutron collision
    Sigma_eth = Xi.*(I_eff + Sigma_sc) ; % (Eq 3.18, Gosse & Phillips, 2001) Effective epithermal loss cross-section (cm2.g-1)
    L_eth = 1./sqrt(3*Sigma_sc.*Sigma_eth); % Epithermal neutron diffusion length (g cm-2)
    L = L_eth;
end
end

