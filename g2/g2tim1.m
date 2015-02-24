function[edt,res] = g2tim1(ww,rhs,grd)
% function [edt,res] = g2tim1(ww,rhs)
% compute edt=1/timestep per cell and abs rel. norm of rhs

% absolute relative residuals
nw = size(ww); res = 0; edt = zeros( nw(2:3) );
for jr=2:nw(3)-1;
  for jp=2:nw(2)-1;
    h1 = abs( rhs(1,jp,jr)./ ww(1,jp,jr) );
    h2 = abs( rhs(2,jp,jr)./(ww(1,jp,jr)+abs(ww(2,jp,jr))) );
    h3 = abs (rhs(3,jp,jr)./(ww(1,jp,jr)+abs(ww(3,jp,jr))) );
    hm = max([h1 h2 h3]);
    res = max(res,hm/grd.dvol(jr)); % relative residual norm
    edt(jp,jr) = hm;
  end
end;
%EOF


