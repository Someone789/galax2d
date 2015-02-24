function dd = g2gradphi(wp,grd,jr,order,doA)
% function dd = g2gradphi(wp,grd,jr,order,doA)
% half gradients for rho,u,v in the phi direction, including
%   gradient w.r.t rho,u,v
% wp is an array at constant R (index jr)
% drp = grd.r(jr) * grd.dphi; % drp is spacing, R*dphi
% if doA: compute Jacobian matrices
% the output dd is a struct

% squared bias for Van Albada limiter
bias = grd.bphi;
if(length(bias) == 1), bias = bias+0*(1:grd.nrad+2); end
% transform wp to {rho,u,v}
er = 1/grd.r(jr);  e1 = 1./wp(1,:);
rh = wp(1,:)*er; u = wp(2,:).*e1; v = wp(3,:).*e1;
% A1 = d {rho,u,v} / d wp
np = length(rh); 
if(doA),
  A1 = zeros(3,3,np); A2a = A1; A2b = A1;
  A1(1,1,:) =     er; 
  A1(2,1,:) = -u.*e1; A1(2,2,:) = e1;
  A1(3,1,:) = -v.*e1; A1(3,3,:) = e1;
else
  A1 = [];
end
%
if(order < 2), % 1st order
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
h2 =   diff(v);
dh = [h1-h2; h1+h2; diff(u)]; % starts at jp=1.5 (between 1 and 2)
% A2a = d(dh)/d ruv_left, A2b = d(dh)/d ruv_right
%  L:  [-4*rhoR/(rhoR+rhoL)^2,0, 1; -4*rhoR/(rhoR+rhoL)^2,0,-1;0,-1,0];
%  R:  [ 4*rhoL/(rhoR+rhoL)^2,0,-1;  4*rhoL/(rhoR+rhoL)^2,0, 1;0, 1,0];
if(doA),
  erq = era.^2; % 4/(rhoR+rhoL)^2
  A2a(1,1,1:end-1) =-erq.*rh(2:end  ); A2a(1,3,1:end-1) =  1;
  A2a(2,1,1:end-1) = A2a(1,1,1:end-1); A2a(2,3,1:end-1) = -1;
    A2a(3,2,1:end-1) = -1;
  A2b(1,1,1:end-1) = erq.*rh(1:end-1); A2b(1,3,1:end-1) = -1;
  A2b(2,1,1:end-1) = A2b(1,1,1:end-1); A2b(2,3,1:end-1) =  1;
    A2b(3,2,1:end-1) =  1;
end
% averaging/limiter towards smallest
dg = zeros(3,np); 
sp = dh.^2+bias(jr); jp=(2:np-1); 
dg(:,jp) = 0.5*(sp(:,jp-1).*dh(:,jp)+...
                sp(:,jp).*dh(:,jp-1))./(sp(:,jp-1)+sp(:,jp));
bm3 = 0.5*(dg(1,jp)+dg(2,jp));
% dbm = d(dg)/d(dh(jp-1)) and dbp = d(dg)/d(dh(jp))
dbm = (dh(:,jp).*dh(:,jp-1)+sp(:,jp  )/2-2*dh(:,jp-1).*dg(:,jp))./(sp(:,jp)+sp(:,jp-1));
dbp = (dh(:,jp).*dh(:,jp-1)+sp(:,jp-1)/2-2*dh(:,jp  ).*dg(:,jp))./(sp(:,jp)+sp(:,jp-1));
% back to rho,u,v
dg = [0.5*rh.*(dg(1,:)+dg(2,:)); dg(3,:); 0.5*(dg(2,:)-dg(1,:))];
%
if(doA),
  BM = zeros(3,3,np); BL = BM; BR = BM; 
  A3 = [1/2 1/2 0;0 0 1;-1/2 1/2 0];
  for jp=2:np-1,
    A3(1,1:2) = rh(jp)/2;
    A3m = A3*diag(dbm(:,jp-1)); 
    A3p = A3*diag(dbp(:,jp-1));
    BL(:,:,jp) = A3m*A2a(:,:,jp-1);
    BR(:,:,jp) = A3p*A2b(:,:,jp  );
    BM(:,:,jp) = A3m*A2b(:,:,jp-1) + A3p*A2a(:,:,jp);
    BM(1,1,jp) = BM(1,1,jp)+bm3(jp-1);
  end
end
% if(0),
%   ih = fix(1+np/2); [BL(:,:,ih) BM(:,:,ih) BR(:,:,ih)]
% end
% periodicity
dg(:,np) = dg(:,   2);
dg(:, 1) = dg(:,np-1);
if(doA),
  BL(:,:,np) = BL(:,:,2); BM(:,:,np) = BM(:,:,2); BR(:,:,np) = BR(:,:,2);
  BL(:,:,1) = BL(:,:,np-1); BM(:,:,1) = BM(:,:,np-1); BR(:,:,1) = BR(:,:,np-1);
else
  BL = []; BM = []; BR = [];
end
% dg: half gradient in rho,u,v;
% BM is d(dg)/d {ruv}M (M, in the centre), similarly for 
% L (left neighbour) and R (right neighbour)
% post-multiply by A1 to get d/dw
dd = struct('rho',rh,'u',u,'v',v,'druv',dg,'A1',A1,'BL',BL,'BM',BM,'BR',BR);
%EOF


