function [qp,xx] = g2plotprepx(ww,grd,aux,iscale,nsmt,xmax,dx,domirror)
% function [qp,xx] = g2plotprepx(ww,grd,aux,iscale,xmax,dx)
% qp = {\log_10 \rho, u, v-v_0} on grid xx,xx
% suggestion: imagesc(xx,xx,qp{1}'); set(gca,'Ydir','normal');
[qq,rha] = g2wtoq(ww,grd,aux,iscale,nsmt);
% log10 density instead of natural log
qq(1,:,:) = qq(1,:,:)/log(10);
if(iscale)
  tsc = ' (scaled per radius)'; % text for title
else
  tsc = ''; % text for title
end
% interpolation
if(nargin > 4), xm = xmax; else xm = gdr.r(end);  end
if(nargin > 5), dd = dx;   else dx = min(grd.dr); end
xm = (-0.5+ceil(xm/dx))*dx; xx = (-xm:dx:xm); yp = (dx/2:dx:xm);
% rh3 = interpn(grd.rm
% variables for plot; remove ghost points in phi
qq = qq(:,2:end-1,:);
%
phi = (0.5+(0:grd.nphi-1))*grd.dphi;
if(nargin > 7), pmirror = domirror; else pmirror = 1; end
[X,Yp] = ndgrid(xx,yp); [P,R]=cart2pol(X,Yp);
qp = cell(3,1);
for k=1:3,
  Fk = griddedInterpolant({phi,grd.r},squeeze(qq(k,:,:)),'linear');
  qk =  Fk(P,R); clear Fk;
  qk = [flipud(fliplr(qk)) qk]; % half space to full space
  if(pmirror), qk = fliplr(qk); end; % mirror
  qp{k} = qk;
  % imagesc(xx,xx,qk'); set(gca,'Ydir','normal');
end
%EOF




