function pbar = g2pbar(gpar,grv,grd)
% function pb = g2pbar(gpar,grv,grd)
% radial dependence of potential V2 for bar
% returns s21p and s31p in struct, where:
%          V2(R,phi) =    V2(R)*cos(2*phi)
%        -dV2/dR     = -s21p(R)*cos(2*phi)
%    -1/R dV2/dphi   = -s31p(R)*sin(2*phi)
%
% gpar.pp     power for spherical density function (typically -1.8)
% gpar.ii     power for windowing function (10 or so)
% gpar.cutoff cut-off if nonzero (1 at co-rotation, 2 at outer Lindblad reso)
pp = gpar.pp; ii = gpar.ii;
% P_2^2(0)*a22*r^(p+2-1), scaled by 1/c^2
vfac22 = 3*grv.a22*grd.r.^(pp+1)/gpar.c^2;
%
nn = 2; n2 = pp+2-nn; n3 = pp+3+nn;
s21p = zeros(size(grd.r)); s31p = s21p;
en23 = 1/(n2*n3);
% note: 3*c22 = grv.b(2) = 3*gpot.a22/(n2*n3);
if(gpar.cutoff == 1 || gpar.cutoff == 2)
  % cut-off bar, density times wm2 := (1-(r/r1)^ii)^2:
  if(gpar.cutoff == 2)
    % cutoff=2: cut-off at the Outer Lindblad Resonance
    % (om0-om)/om1 = -1/2, om0 = fa0*r^(p/2), om1 = fa1*om0
    fa1 = sqrt(pp+4); r1 = ( gpar.om/(grv.fa0*(1+0.5*fa1)) )^(2/pp);
  else
    % cutoff=1: cut-off at co-rotation
    r1 = grv.rcr;
  end
  zeta = grd.r/r1;
  % r > r1:
  %    2*anm*ii^2*zeta^(-n3)/((2*n+1)*(p+n+3)*(p+n+3+ii)*(p+n+3+2*ii))
  %     *3*grd.r.^(pp+1)/gpar.c^2
  %    only consider n=m=2
  jout = find(zeta >= 1);
  if(~isempty(jout)),
    % note: r^(p+1-n3)=r^(-4)
    c22o = 2*ii^2/((2*nn+1)*n3*(n3+ii)*(n3+2*ii));
    s31p(jout) = (2*c22o)*vfac22(jout).*zeta(jout).^(-n3);
    s21p(jout) = 1.5*s31p(jout);
  end
  % r < r1:
  jin = setdiff((1:length(grd.r)),jout);
  if(~isempty(jin)),
    zin = zeta(jin); ct = zeros(4,length(zin));
    ct(1,:) =  1/( n2*       n3); % = en23
    ct(2,:) = -2/((n2+  ii)*(n3+  ii)) * zin.^ii;
    ct(3,:) =  1/((n2+2*ii)*(n3+2*ii)) * zin.^(2*ii);
    ct(4,:) = -2*ii^2/((2*nn+1)*n2*(n2+ii)*(n2+2*ii))*zin.^(-n2);
    sumct = sum(ct,1);
    s31p(jin) = -2*sumct.*vfac22(jin);
    cu = ct; fc = [pp+2,pp+2+ii,pp+2+2*ii,2];
    for k=1:4, cu(k,:) = ct(k,:)*fc(k); end
    sumcu = sum(cu,1);
    s21p(jin) =   sumcu.*vfac22(jin);
  end
else % no cut-off;same as r<r1 with ct(2:4,:)=0 above
  % c_{nm}*r^(p+2)*cos(m*phi), c_{nm} = a_{nm}/((p+3+n)*(p+2-n)):
  s21p =  (pp+2)*en23*vfac22;
  s31p =      -2*en23*vfac22;
  r1 = 0;
end
% save into struct
%       V2(R,phi)=    V2(R)*cos(2*phi)
%     -dV2/dR    = -s21p(R)*cos(2*phi)
% -1/R dV2/dphi  = -s31p(R)*sin(2*phi)
pbar = struct('s21p',s21p,'s31p',s31p,'r1',r1);
%EOF

