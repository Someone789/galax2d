function [qq,rha] = g2wtoq(ww,grd,aux,iscale,nsmt)
% [qq,rhoa] = g2wtoq(ww,grd,aux,iscale)
% converts ww = R*rho*[1,u,v] to qq = {log(rho),u,v-v0}
% ghost cells are not removed
% if iscale > 0, rhoa is average density over phi per radius R
%   optionally smoothed nmst times by convolution with [1/4,1/2,1/4]
qq = ww; qq(2,:,:) = qq(2,:,:)./qq(1,:,:); qq(3,:,:) = qq(3,:,:)./qq(1,:,:);
nr = size(qq,3);
for jr=1:nr,
  qq(1,:,jr) = qq(1,:,jr)./grd.r(jr);  % times 1/R
  qq(3,:,jr) = qq(3,:,jr)-aux.v0m(jr); % subtract circular velocity
end
if(nargin > 3 & iscale), 
  % divide by average density over phi
  for jr=1:nr, rha(jr) = mean(qq(1,2:end-1,jr)); end;
  if(nargin > 4 & nsmt), % some smoothing?
    ir = (1:nr);
    for k=1:nsmt,
      rha = 0.5*(rha(ir)+0.5*(rha(min(nr,ir+1))+rha(max(1,ir-1))));
    end
  end
  for jr=1:nr, qq(1,:,jr) = qq(1,:,jr)/rha(jr); end
else
  rha = [];
end
% log density (not log10!)
qq(1,:,:) = log(qq(1,:,:));
%EOF
