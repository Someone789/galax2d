function grd = g2bias(ww,grd,aux)
% function grd = g2bias(ww,grd,aux)
% determines squared bias in radial and phi direction for limiter
bphi = zeros(1,grd.nrad+2);
% _____________________________________________________________________________
% phi-direction
er = 1./grd.r;
for jr=2:grd.nrad+1,
  wp = ww(:,1:end-1,jr); e1 = 1./wp(1,:); 
  rh = wp(1,:)*er(jr); u = wp(2,:).*e1; v = wp(3,:).*e1;
  era = 2./(rh(1:end-1)+rh(2:end));
  h1 = diff(rh).*era; h2 = diff(v);
  dh = [h1-h2; h1+h2; diff(u)]; 
  bphi(jr) = sum(sum(abs(dh)))/(3*grd.nphi);
end  
% _____________________________________________________________________________
% R-direction
da = zeros(3,grd.nrad+2);
for jp=2:grd.nphi+1,
  wp = ww(:,jp,:); % wp(:,1) = wp(:,2); wp(:,end-1) = wp(:,end); 
  e1 = 1./wp(1,:);
  rh = wp(1,:).*er; u = wp(2,:).*e1; vv = wp(3,:).*e1; v = vv-aux.v0m;
  era = 2./(rh(1:end-1)+rh(2:end));
  h1 =  diff(rh).*era; h2 = diff(u);
  dh = abs( [h1-h2;  h1+h2; diff(v)] );
  for j=2:grd.nrad+1, 
    for k=1:3, da(k,j) = da(k,j)+min(grd.hrl(j)*dh(k,j-1),grd.hrr(j)*dh(k,j)); end
  end
end
brad = sum(da)/(3*grd.nrad);
% smoothing
if(1),
  brad = smooth1(brad);
  bphi = smooth1(bphi);
end
% now square the result
grd.brad = brad.^2; grd.bphi = bphi.^2;
% _____________________________________________________________________________
function a = smooth1(a)
% function a = smooth1(a)
% simple smoothing by convolution with {1/4,1/2,1/4}
na = length(a); jm = [1 (1:na-1)]; jp = [(2:na) na];
a = 0.5*(a+0.5*(a(jp)+a(jm)));
%EOF
