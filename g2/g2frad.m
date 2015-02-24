function [f,Af] = g2frad(isg,R,dphi,rho,uu,vv)
% function [f,Af] = g2frad(isg,R,dphi,rho,uu,vv);
% flux in R-direction, f = dphi*R*rho*[u;u^2+c^2;u*v]; scaled to c=1
%   isg: +1 or -1 for + or - flux
% Af is the Jacobian w.r.t. {rho,u,v}
% dphi is the grid spacing in the phi-direction
us = isg*uu;
if(us <= -1),
  f = zeros(3,1); Af = zeros(3,3);
else
  if(us >= 1)
    % full flux and jacobian w.r.t. rho,uu,vv; c is scaled to 1
    w1b = dphi*R*rho; f = w1b*[uu; uu*uu+1; uu*vv];
    Af = [uu     ,      rho,      0; ...
          uu*uu+1, 2*rho*uu,      0; ...
          uu*vv  ,   rho*vv, rho*uu];
  else 
    % |us| < 1, split flux:
    % f = R*r/(4*cs)*(u+cs)^2*[1;2*cs;v], cs = isg*c, c=1
    w1b = dphi*R*rho; uc = uu+isg;
    is2 = 0.5*isg; is4 = 0.25*isg;
    f = w1b*is4*uc^2*[1; 2*isg; vv];
    Af = (is2*uc)*[0.5*uc,         rho, 0; ...
                   isg*uc,    2*isg*rho, 0; ...
                   0.5*uc*vv,    vv*rho, 0.5*uc*rho];
  end
  Af = (dphi*R)*Af;
end
%EOF
