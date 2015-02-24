function mm = g2mesh(pp,kappa, Rmin,Rmax, nrad,nphi)
% function mm = g2mesh(pp,kappa, Rmin,Rmax, nrad,nphi)
% create mesh struct "mm"
% nrad cells in radial directions between Rmin and Rmax, 
%   2 ghost cells will be added at each side
% nphi cells in angular direction phi
% equidistant in xi=R^(1+kappa*pp/2)
ph1 = 1+kappa*pp*0.5; xi0 = Rmin^ph1; xi1 = Rmax^ph1; dxi = (xi1-xi0)/nrad;
rl = (xi0+(0:nrad)*dxi).^(1/ph1); % left side of cell
rl(1) = Rmin; rl(end) = Rmax; % round-off correction
rl = [Rmin rl Rmax]; % add dummy entries
r = 0.5*(rl(1:end-1)+rl(2:end)); % 1/2*(rl+rr)
dr = diff(rl); dr(1) = 0; dr(end) = 0; 
dphi = (pi/nphi);
dvol = dr*dphi; % dR dphi
% squared bias for limiter
brad = (0.01*dr).^2;
bphi = zeros(size(brad))+(0.01*dphi).^2; 
%
mm = struct('rl',rl,'r',r,'dr',dr,'dphi',dphi,'dvol',dvol,...
            'nrad',nrad,'nphi',nphi,'hrl',[],'hrr',[],...
            'bphi',bphi,'brad',brad);
mm = g2meshaver(mm);
%EOF