function grc = g2meshcoarsen(grd,nc)
% subsample mesh in grd to nc cells
grc = grd; nrat = grd.nrad/nc;
if(nrat <= 1), return; end
grc.dphi = nrat*grd.dphi;
% new cell boundaries in R
grc.rl = [grd.rl(1) grd.rl(2:nrat:end-1) grd.rl(end)];
grc.r  = 0.5*(grc.rl(1:end-1)+grc.rl(2:end));
grc.dr = diff(grc.rl); grc.dr(1) = 0;  grc.dr(end) = 0; 
grc.dvol = grc.dr*grc.dphi; % dR dphi
grc.nrad = nc; grc.nphi = nc; 
grc.brad = (0.01*grc.dr).^2;
grc.bphi = zeros(size(grc.dr))+(0.01*grc.dphi).^2; 
%
grc = g2meshaver(grc);
%EOF

