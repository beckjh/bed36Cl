function sc = scaling(perp_below_coll,heights,thick,Lambda_c,Lambda_r,alpha,beta,gamma,rho_colls,rho_rock)
% Calculate the scaling factor sc associated with directions passing the
% colluvial wedge.
%   perp_below_colll: perpendicular distance of sample to colluvial surface
%   heights: distance to gamma dip along beta plane (cm)
assert(beta>0&& beta<90);assert(alpha>=0 && alpha<beta);
d_thick=0.1;
d_theta = 0.1;
d_phi = 0.1;
thicks = d_thick:d_thick:thick;
thicks = reshape(thicks,1,1,[]);
depths=perp_below_coll/sind(beta-alpha);%along beta plane
assert(all(depths>=0));
m = 2.3 ;
alpha = alpha*pi/180 ;
beta = beta*pi/180 ;
gamma = gamma*pi/180;
thetas = (0:d_theta:pi/2)';
phis = 0:d_phi:2*pi;
sc=zeros(numel(depths),1);
for i = 1:numel(depths)
    depths_diffs = bsxfun(@rdivide,thicks,tan(pi/2-beta+thetas));% N_1 x 1 x N_3
    depths_under_colluvium=max(0,bsxfun(@minus,depths(i),depths_diffs));%N_1 x 1 x N3
    perp1 = depths_under_colluvium*sin(beta-alpha);%perpendicular length through colluvial wedge, N_1 x 1 x N_3 (same as input, but corrected for thickness)
    d1 = max(0,bsxfun(@times,f(alpha,thetas,phis),perp1));%N_1 x N_2 x N_3
    through_colluvial = g(beta,thetas,phis,depths(i),thicks);
    d1(~through_colluvial)=0;
    perp2 = thicks;%perp. distance to beta part of rock 1 x 1 x N_3
    no_gamma=g(beta,thetas,phis,heights(i),thicks);
    tt=f(beta,thetas,phis);
    tt(tt<0)=Inf;
    d2_nogamma = bsxfun(@times,tt,perp2);% N_1 x N_2 x N_3 travel distance to beta part
    h_below_gamma_dip = bsxfun(@plus,heights(i)*sin(beta),thicks*cos(beta));%1 x 1 x N3 
    l_in_front_gamma_dip = bsxfun(@minus,heights(i)*cos(beta),thicks*sin(beta));%1 x 1 x N3 
    perp3 = cos(gamma)*(h_below_gamma_dip - tan(gamma)*l_in_front_gamma_dip);% 1 x 1 x N3 perp. distance to gamma part of rock
    tt=f(gamma,thetas,phis);
    tt(tt<0)=Inf;
    d2_gamma = bsxfun(@times,tt,perp3);%full size
    d2=zeros(size(d2_gamma));
    d2(no_gamma)=d2_nogamma(no_gamma);
    d2(~no_gamma)=d2_gamma(~no_gamma);
    attenuation = bsxfun(@times,exp(-rho_colls(i)*d1./Lambda_c(i)),exp(-rho_rock*d2./Lambda_r(i)));%full size
    I = bsxfun(@times,cos(thetas).^m,attenuation);%full size
    dA = abs(sin(thetas))*d_phi*d_theta/numel(thicks);% N_1
    sc(i) = sum(sum(sum(bsxfun(@times,I,dA),1),2),3)*(m + 1)/(2*pi) ;% 1 x 1 x 1
end
end

function d = f(dip,theta,phi)%N1 x N2 x 1 x 1
d=1./bsxfun(@plus,cos(dip)*cos(theta),bsxfun(@times,sin(phi)*sin(dip),sin(theta)));%abs only necessary due to roundoff. (Affects only very large lengths, which will then be practically ignored)
end

function hit = g(dip,theta,phi,d,thick)
    h = bsxfun(@plus,d*sin(dip),thick*cos(dip));%1 x 1 x N3 
    l = bsxfun(@minus,d*cos(dip),thick*sin(dip));%1 x 1 x N3 
    c = bsxfun(@rdivide,l,abs(sin(phi)));
    h_cross = bsxfun(@rdivide,c,tan(theta));
    hit=NaN*zeros(numel(theta),numel(phi),numel(thick));
    hit(:,phi<=pi,l>=0) =true;
    hit(:,phi>pi,l<0) = false;
    h_cross_larger_than_h=bsxfun(@gt,h_cross,h);
    hit(:,phi>pi,l>=0) =  h_cross_larger_than_h(:,phi>pi,l>=0);
    hit(:,phi<=pi,l<0) = ~h_cross_larger_than_h(:,phi<=pi,l<0);
    hit=logical(hit);
end
