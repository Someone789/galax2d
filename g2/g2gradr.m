function dd = g2gradr(wp,grd,aux,order,doA)
% function dd = g2gradr(wp,grd,aux,order,doA)
% half gradients for rho,u,v in the R direction, including
%  gradient w.r.t rho,u,v
% wp is an array at constant constant phi (index jp) with 2 ghost cells
% requires aux.v0m (axisymmetric solution for v)
% if doA: compute Jacobian matrices
% the output dd is a struct

% squared bias for Van Albada limiter
bias = grd.brad;
% transform wp to {rho,u,v}
er = 1./grd.r;  e1 = 1./wp(1,:);
rh = wp(1,:).*er; u = wp(2,:).*e1; vv = wp(3,:).*e1;
v = vv-aux.v0m; % tangential velocity relative to circular case
% A1 = d {rho,u,v} / d wp
np = length(rh);
if(doA),
  A1(1,1,:) =      er;
  A1(2,1,:) =  -u.*e1; A1(2,2,:) = e1;
  A1(3,1,:) = -vv.*e1; A1(3,3,:) = e1;
else
  A1 = [];
end
% test 1st order
if(order < 2),
  if(doA), BM = zeros(3,3,np); else BM = []; end
  BL = BM; BR = BM;
  dg = zeros(3,np);
  dd = struct('rho',rh,'u',u,'v',v,'druv',dg,...
              'A1',A1,'BL',BL,'BM',BM,'BR',BR);
  return;
end
% transform to characteristic variables
era = 2./(rh(1:end-1)+rh(2:end));
h1 =   diff(rh).*era;
h2 =   diff(u);
dh = [h1-h2;  h1+h2; diff(v)]; % starts at jr=1.5 (between 1 and 2)
% _____________________________________________________________________________
% boundaries: extrapolate from interior if outflow occurs
iexl = [0 0 0]; iexr = [0 0 0]; % fill if extrapolated
%
doextrapolate = 0;
if(doextrapolate == 0),
  hr1 = 0; hrn = 0; % do not extrapolate
else
  % extrapolation based on characteric variables
  % causes convergence problems on 512^2; bugs?
  hr1 = grd.hrr(1); hrn = grd.hrl(np);
  if(u(   2)  <  1), % inner boundary
    dh(1,   1) = hr1*dh(1,   2); iexl(1) = 1; % -char
    if(u(  2) <  0), dh(3,   1) = hr1*dh(3,   2); iexl(3) = 1; end % 0char
    if(u(  2) < -1), dh(2,   1) = hr1*dh(2,   2); iexl(2) = 1; end % +char
  end
  if(u(np-1) >= -1), % outer boundary
    dh(2,np-1) = hrn*dh(2,np-2); iexr(2) = 1; % +char
    if(u(np-1)>= 0), dh(3,np-1) = hrn*dh(3,np-2); iexr(3) = 1; end % 0char
    if(u(np-1)>= 1), dh(1,np-1) = hrn*dh(1,np-2); iexr(1) = 1; end % -char
  end
end
% _____________________________________________________________________________
% A2a = d(dh)/d ruv_left, A2b = d(dh)/d ruv_right
% L:  [-4*rhoR/(rhoR+rhoL)^2, 1,0; -4*rhoR/(rhoR+rhoL)^2,-1,0; 0 0 -1];
% R:  [ 4*rhoL/(rhoR+rhoL)^2,-1,0;  4*rhoL/(rhoR+rhoL)^2, 1,0; 0 0  1];
% (1)|(2)...(np-1)|(np)  -- rho
%    1   2       np-1    -- A2a,b
if(doA),
  erq = era.^2; % 4/(rhoR+rhoL)^2
  A2a = zeros(3,3,np-1); A2b = A2a;
  A2a(1,1,:) =-erq.*rh(2:end); A2a(1,2,:) =  1;
  A2a(2,1,:) = A2a(1,1,:); A2a(2,2,:) = -1;
  A2a(3,3,:) = -1;
  A2b(1,1,:) = erq.*rh(1:end-1); A2b(1,2,:) = -1;
  A2b(2,1,:) = A2b(1,1,:); A2b(2,2,:) =  1;
  A2b(3,3,:) =  1;
end
% averaging/limiter towards smallest
dg = zeros(3,np);
dbm = zeros(3,np-2); dbp = dbm; bm3 = zeros(1,np-2);
for jp=2:np-1
  h1l = grd.hrl(jp)*dh(:,jp-1); sm = h1l.^2+bias(jp);
  h1r = grd.hrr(jp)*dh(:,jp  ); sp = h1r.^2+bias(jp);
  dg(:,jp) = 0.5*(sm.*h1r+sp.*h1l)./(sm+sp);
  % dbm = d(dg)/d(dh(jp-1)) and dbp = d(dg)/d(dh(jp))
  dbm(:,jp-1) = grd.hrl(jp).*(sp/2+h1l.*(h1r-2*dg(:,jp)))./(sm+sp);
  dbp(:,jp-1) = grd.hrr(jp).*(sm/2+h1r.*(h1l-2*dg(:,jp)))./(sm+sp);
  % d 0.5*rh(jp)*(da(1)+da(2))/d rho
  bm3(jp-1) = 0.5*(dg(1,jp)+dg(2,jp));
end
%
for jp=2:np-1
  % back to rho,u,v
  da = dg(:,jp);
  dg(:,jp) = [0.5*rh(jp)*(da(1)+da(2)); 0.5*(da(2)-da(1)); da(3)];
end
%
if(doA),
  BM = zeros(3,3,np); BL = BM; BR = BM;
  % endpoints:
  A2a(:,:,   1) = 0;
  A2b(:,:,np-1) = 0;
  % A3 handles step from da to dg
  A3 = [1/2 1/2 0;-1/2 1/2 0;0 0 1];
  for jp=2:np-1,
    A3(1,1:2) = rh(jp)/2;
    A3m = A3*diag(dbm(:,jp-1));
    A3p = A3*diag(dbp(:,jp-1));
    BL(:,:,jp) = A3m*A2a(:,:,jp-1);
    BM(:,:,jp) = A3m*A2b(:,:,jp-1) + A3p*A2a(:,:,jp);
    BR(:,:,jp) =                     A3p*A2b(:,:,jp);
    % extrapolate
    if(doextrapolate),
      if(jp ==    2), % first interior cell: A3p=A3p+hr1*A3m;A3m=0; if iexl(k)
        % BL(:,:,jp) = 0; % already zero via A2a(:,:,1)
        % P2b = A2b(:,:,jp-1), expected to be 0
        P2b = A2b(:,:,jp-1);
        Q2a = zeros(3,3); Q2b = zeros(3,3);
        for k=1:3,
          if(iexl(k)),
            P2b(k,:) = 0; Q2a(k,:) = A2a(k,:,jp); Q2b(k,:) =  A2b(k,:,jp);
          end
        end
        BM(:,:,jp) = A3m*(P2b+hr1*Q2a) + A3p*A2a(:,:,jp);
        BR(:,:,jp) = A3m*(    hr1*Q2b) + A3p*A2b(:,:,jp);
      end
      %
      if(jp == np-1), % last interior cell: A3m=A3m+hrn*A3p;A3p=0; if iexl(k)
        % BR(:,:,jp) = 0; % already zero via A2b(:,:,np-1)
        P2a = A2a(:,:,jp); % P2a = A2a(:,:,jp), expected to be 0
        Q2a = zeros(3,3);  Q2b = zeros(3,3);
        for k=1:3,
          if(iexr(k)),
            P2a(k,:) = 0; Q2a(k,:) = A2a(k,:,jp-1); Q2b(k,:) = A2b(k,:,jp-1);
          end
        end
        BL(:,:,jp) = A3m*A2a(:,:,jp-1) + A3p*(    hrn*Q2a);
        BM(:,:,jp) = A3m*A2b(:,:,jp-1) + A3p*(P2a+hrn*Q2b);
      end
    end
    % contribution from d [ 0.5*rh(jp)*(da(1)+da(2)) ]/d rho
    BM(1,1,jp) = BM(1,1,jp)+bm3(jp-1);
  end
else
  BL = []; BM = []; BR = [];
end
% dg: half gradient in rho,u,v;
% BM is d(dg)/d {ruv}M (M, in the centre), similarly for
% L (left neighbour) and R (right neighbour)
% please post-multiply by A1 to get d(.)/dw
dd = struct('rho',rh,'u',u,'v',v,'druv',dg,...
  'A1',A1,'BL',BL,'BM',BM,'BR',BR);
%EOF
