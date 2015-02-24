function [g,Ag] = g2fphi(isg,dR,rho,uu,vv)
% function [g,Ag] = g2fphi(isg,dR,rho,uu,vv);
% flux in phi-direction, g = rho*[v;uv;v^2+c^2]; scaled to c=1
%   isg: +1 or -1 for + or - flux
% Ag is the Jacobian w.r.t. {rho,u,v}
% dR is the grid spacing in the radial R-direction
vs = isg*vv; 
if(vs <= -1),
  g = zeros(3,1); Ag = zeros(3,3);
else
  if(vs >= 1)
    % full flux and jacobian w.r.t. rho,uu,vv; c is scaled to 1 
    g  = (dR*rho)*[vv; vv*uu; vv*vv+1];
    Ag = [vv    ,      0,      rho; ...
          uu*vv , rho*vv,   rho*uu; ...
          vv^2+1,      0, 2*rho*vv];
  else 
    % if(vs > -1), 
    % split flux: dR*g = dR*rho/(4*cs)*(v+cs)^2*[1;u;2*cs], cs = isg*c, c=1
    vc = vv+isg; is2 = 0.5*isg; is4 = 0.25*isg; 
    g  = (dR*rho*is4*vc^2)*[1; uu; 2*isg];
    Ag = is2*vc*[0.5*vc   , 0         ,       rho; ...
                 0.5*vc*uu, 0.5*vc*rho,    uu*rho; ...
                 isg*vc   , 0         , 2*isg*rho];
  end
  Ag = dR*Ag;  % note that 1/R sits in d{ruv}/dw
end


