function ppm=ppm_function(sample,A_k)
% Conversion of oxyde percents into percents of the oxyded element
chimie = sample(:,1:62) ;
ppm = chimie ;
ppm(:,44) = chimie(:,44)*A_k(44)/(A_k(44) + 2*A_k(59)) ; % Si in percent
ppm(:,45) = chimie(:,45)*2*A_k(45)/(2*A_k(45) + 3*A_k(59)) ; % Al in percent
ppm(:,46) = chimie(:,46)*2*A_k(46)/(2*A_k(46) + 3*A_k(59)) ; % Fe in percent
ppm(:,47) = chimie(:,47)*A_k(47)/(A_k(47) + A_k(59)) ; % Mn in percent
ppm(:,48) = chimie(:,48)*A_k(48)/(A_k(48) + A_k(59)) ; % Mg in percent
ppm(:,49) = chimie(:,49)*A_k(49)/(A_k(49) + A_k(59)) ; % Ca in percent
ppm(:,50) = chimie(:,50)*2*A_k(50)/(2*A_k(50) + A_k(59)) ; % Na in percent
ppm(:,51) = chimie(:,51)*2*A_k(51)/(2*A_k(51) + A_k(59)) ; % K in percent
ppm(:,52) = chimie(:,52)*A_k(52)/(A_k(52) + 2*A_k(59)) ; % Ti in percent
ppm(:,53) = chimie(:,53)*2*A_k(53)/(2*A_k(53) + 5*A_k(59)) ; % P in percent
ppm(:,56) = chimie(:,56)*2*A_k(56)/(2*A_k(56) + A_k(59)) ; % H water in percent
O_water = chimie(:,56) - ppm(:,56) ; % O_water in percent
ppm(:,58) = chimie(:,58)*A_k(58)/(A_k(58) + 2*A_k(59)) ; % C in percent
ppm(:,59) = sum([chimie(:,44:53) chimie(:,58)],2) - sum([ppm(:,44:53) ppm(:,58)],2) ; % O rock in percent
ppm(:,60) = O_water ;
ppm(:,44:53) = ppm(:,44:53)*1e+4 ; % in ppm
ppm(:,56) = ppm(:,56)*1e+4 ; % in ppm
ppm(:,58) = ppm(:,58)*1e+4 ; % in ppm
ppm(:,59) = ppm(:,59)*1e+4 ; % in ppm
ppm(:,60) = ppm(:,60)*1e+4 ; % in ppm
end