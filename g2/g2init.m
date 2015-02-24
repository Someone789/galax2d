function [w,aux] = g2init(gpar,grv,grd)
% function [w,aux] = g2init(gpar,grv,grd)
% Initialize solution vector w in each cell.
% struct aux contains auxiliary arrays:
%   Axisymmetric velocity v0l at left (or lower) side of cell and
%   v0m as average per cell, based on R*v0 with volume element R dR.
%   Central solution cen = [rho,u,v] as a function of phi

% _____________________________________________________________________________
% density: initial, inner boundary, outer boundary
rhoinit = gpar.rhoinit; rhoinner = gpar.rhoinner; rhoouter  = gpar.rhoouter;
% _____________________________________________________________________________
% rotation: vrot = v+om*R = om0*R = fa0*R^(1+p/2)
%           v0/c = vrot/c-om/c*R = R*(  (fa0/c)*R^(p/2) - (om/c) )
%           with fa0 = sqrt(c00*(p+2));
% integrate R*v/c w.r.t. R to get v0m
  omc = gpar.om/gpar.c; % omega/c
  sa0 = grv.fa0/gpar.c; % fa0/c
  pp  = gpar.pp; pph = 0.5*pp; % pph=p/2
dr = diff(grd.rl); % delta R
% v0l is velocity at left (or lower) side of cell,
% v0m is the average over the cell, using R*v0 with volume element R dR
v0l = grd.rl   .*( sa0        *grd.rl.^pph -  omc   );
v0m = grd.rl.^3.*( sa0/(3+pph)*grd.rl.^pph -  omc/3 ); v0m = diff(v0m);
% divide by R*dR, except at endpoints where dr=0
v0m(2:end-1) = v0m(2:end-1)./(grd.r(2:end-1) .* dr(2:end-1));
% endpoints
v0m(1) = v0l(2); v0m(end) = v0l(end-1); % endpoints
% _____________________________________________________________________________
% fill conserved quantities: w = R*rho*{1,u,v}
% note: one extra cell at each side, whereas the mesh has two ...
w = zeros(3,grd.nphi+2,grd.nrad+2);
for j=2:grd.nrad+1,
  w(1,:,j) = rhoinit*grd.r(j);
  w(3,:,j) = w(1,:,j)*v0m(j);
end
% _____________________________________________________________________________
% log(rho/rhoc) = eps*R^(  p/2) A0 cos(2*phi)
%  u            = eps*R^(1+p/2) A1 sin(2*phi)
%  v-v0m        = eps*R^(1+p/2) A2 cos(2*phi)
% note: total v requires extra v0R = (om0-omc)*R
% om0 = sqrt(V_R/R) = vrot/R = sqrt( (p+2)*c0*R^p )
ep0 = grv.bc(2)/sa0; % eps = c2(R)/sqrt((p+2)*c0)
r0  = grd.r(1); % inner radius
om0 = sa0*r0^pph;
om1 = sqrt(pp+4)*om0; h0 = om0-omc; h2 = om1^2-4*h0^2;
% amplitudes when cq/R^2 is neglected
A2 = ( (1+pph)*om1^2 + 4*h0*om0 )/h2;
A1 = om0/h0*( A2-(1+pph) );
A0 = pp*om0/h0*( (om0-2*omc)/h2 + 8*omc*h0*((1+pph)*h0+om0)/h2^2 );
%
aa1 = ep0*A0*r0^pph;
aa2 = ep0*A1*r0^(1+pph);
aa3 = ep0*A2*r0^(1+pph);
%
cen = zeros(3,size(w,2)); % array in phi direction
tpi = (-1:2:2*grd.nphi+1)*grd.dphi; % value of 2*phi inside cells, 2 extra
cosp = cos(tpi); sinp = sin(tpi);
cen(1,:) = r0*exp(aa1*cosp); % to be multiplied by rhoinner
cen(2,:) =        aa2*sinp;  % u1
cen(3,:) = v0m(1)+aa3*cosp;  % v0+v1
% w1=rhoinner*cen1;w2=w1*cen2;w3=w1*cen3 for jr=1, see g2wfix
% _____________________________________________________________________________
% sin(2*phi) and cos(2*phi)
dphi = pi/grd.nphi; % phi = (-dphi:2*dphi:2*pi+dphi);
phi = (-1:2:2*grd.nphi+1)*dphi;
sinp = sin(phi); cosp = cos(phi);
% _____________________________________________________________________________
% V = R^(p+2)*c0 + ep*rc2(R)*cos(2*phi); rc2(R) = c2*R^(p+2) or c2(R)*R^(p+2)
% s = rho {0, -R dV/dR + c^2 + (v+om R)^2, -dV/d\phi-u(v+2 om R)}
% ds/dw=[0,0,0;
%    -diff(V,R)+om^2*R+c^2/R-vv^2/R, 0, 2*(om+vv/R);
%    uu*vv/R-diff(V,phi)/R, -(2*om+vv/R), -uu/R]
%
% (2,1): V = V0+O(ep), -R*diff(V0,R) = -(om0*R)^2
%        v0 = (om0-om)*R, v = v0+O(ep),
%   sw21 = b21-vv^2/R, b21 = c^2+(om^2-om0^2)*R - ep*rc2'*cos(2*phi)
%  -diff(V,R)+om^2*R = (om^2-om0^2)*R = 
%                    = -R*(om0-om)*(om0+om) = -v0*(v0+2*om*R)/R
% bhh*c^2 = 1/R-v0*(v0+2*om*R)/R
bhh = (1-v0m.*(v0m+2*omc*grd.r))./grd.r;
% radial dependence of potential of bar
pbar = g2pbar(gpar,grv,grd);
% contribution to jacobian of source term, ds/dw
% save storage by not using b21 and b31?
b21 = zeros(grd.nrad+2,grd.nphi+2); b31 = b21;
for jr=2:grd.nrad+1,
  for jp=2:grd.nphi+1,
     b21(jp,jr) = bhh(jr)-pbar.s21p(jr)*cosp(jp);
     b31(jp,jr) =        -pbar.s31p(jr)*sinp(jp);
  end;
end;
% _____________________________________________________________________________
aux = struct('v0m',v0m,'v0l',v0l,'cen',cen,'sinp',sinp,'cosp',cosp,...
             'pbar',pbar,'bhh',bhh,'b21',b21,'b31',b31);
% set boundary values
w = g2wfix(w,aux,gpar,grd);
%EOF




