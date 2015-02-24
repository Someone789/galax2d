function grd = g2meshaver(grd)
% weights for averaging/limiter on non-equidistant grid
hrl = zeros(size(grd.r)); hrr = hrl;
% hl(j) = 2*d(j)/(d(j)+d(j-1)), hr(j) = 2*d(j)/(d(j)+d(j+1))
% (d(j)+d(j-1))/2 = r(j)-r(j-1)
  drq = (1/12)*(grd.dr./grd.r).^2; % 3rd order correction
  drq = 0*drq;
  hl = grd.dr+2*grd.r.*drq;
  hr = grd.dr-2*grd.r.*drq;
  hd = 1./diff(grd.r.*(1+drq));
  hrl(2:end-1) = hl(2:end-1).*hd(1:end-1);
  hrr(2:end-1) = hr(2:end-1).*hd(2:end);
  % boundaries: factors for extrapolation
  hrl(  1)=0; hrr(  1)=0.5*hl(    2)*hd(    2);
  hrr(end)=0; hrl(end)=0.5*hr(end-1)*hd(end-1);
%
grd.hrr = hrr; grd.hrl = hrl;
%EOF
