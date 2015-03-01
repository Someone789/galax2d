function cm3 = g2color(type,nshift,power)
% function cm3 = g2color(type,nshift,power)
% define a colormap, type = 0 is gray, 1 is color
% defaults: 1,6,1
% or use colormapeditor ...
if(nargin < 1), type = 1; end
if(type == 0),
  if(nargin < 2), aw = 2; else aw = nshift; end
  h0 = (127:-1:0)'/127; h0 = 2*h0-1;
  if(aw == 0), h1 = h0; 
  elseif(aw < 0), h1 = atan(h0*abs(aw));  
  else h1 = tanh(aw*h0); 
  end
  h1 = h1/max(abs(h1)); h1 = 0.5*(1+h1);
  cm3 = [h1,h1,h1]; colormap(cm3);
  return;    
end
if(nargin > 1), nshft = nshift; else nshft = 6; end
if(nargin > 2), gpow = power; else gpow = 1; end
% half the number of colors
nlen = 64;
%
nle1 = nlen-1; rm0 = (0:nle1)/nle1;
rm1 = [(0:nle1-1)/(nle1+nshft) 1]; 
rm3 = max(0,min(1,rm1)); 
if((gpow ~= 1) && (gpow ~= 0)), rm3 = rm3.^gpow; end
% negative values:
r1 = rm1.^2; 
g1 = r1; b1 = max(r1,min(1,6*r1));
% positive red,green.blue
r2 = 1+0*rm3; rmf = fliplr(rm3); g2 = rmf; b2 = rmf;
g2 = g2-0.1+0.65*sin(pi*rm0).^4; 
b2 = b2+2*sin(pi*max(0,rm0-0.7)).^4;
% clip
g2 = max(rmf,min(1,g2));
b2 = max(rmf,min(1,b2));
%
cm3 = [r1 r2; g1 g2; b1 b2]'; cm3 = max(0,min(1,cm3)); colormap(cm3);
%EOF
