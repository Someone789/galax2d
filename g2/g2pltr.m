function g2pltr(ww,grd,aux,conv,iscale)
% function g2pltr(ww,grd,aux,conv,iscale)
% plot ww, scaled by average per ring (if iscale=1)
% ww(1,:,:) = rho*R as function of phi and R
% ww(2,:,;) = w1*u, ww(3,:,:) = w1*v
% conv is convergence history

persistent myFigureHandle;
if(isempty(myFigureHandle) || ~ishandle(myFigureHandle)),
  myFigureHandle = gcf;
end

rh = squeeze(ww(1,:,:)); nr = size(rh,2); rha = zeros(nr,1);
% mean of rho per radius 
for jr=1:nr, rha(jr) = mean(rh(2:end-1,jr)); end;
%
if(nargin > 4 && iscale == 1),
  if(0), % some smoothing?
    ir = (1:nr);
    rha = 0.5*(rha(ir)+0.5*(rha(min(nr,ir+1))+rha(max(1,ir-1))));
  end
  % divide (note that R occurs in numerator and denominator)
  for jr=1:nr, rh(:,jr) = rh(:,jr)/rha(jr); end;
  tsc = ' (scaled per radius)'; % text for title
else
  % divide by mean R=(Rl+Rr)/2
  for jr=1:nr, rh(:,jr) = rh(:,jr)/grd.r(jr); end;
  tsc = ''; % text for title
end
% mean density
rha = rha./grd.r';
%
phi = (0:2*grd.nphi)*grd.dphi; rr = grd.rl(2:end-1);
% patches
nphi = length(phi); nrr = length(rr); vrt = zeros(nphi,nrr,2);
for jr=1:nrr,
  for jp=1:nphi,
    vrt(jp,jr,:) = rr(jr)*[cos(phi(jp)) sin(phi(jp))];
  end;
end
pat = zeros((nphi-1),(nrr-1),4);
for jr=1:nrr-1,
  for jp=1:nphi-1, kk = jp+nphi*(jr-1) ;
    pat(jp,jr,:) = [kk kk+1 kk+1+nphi kk+nphi];
  end;
end
% variables for plot; remove ghost cells
rh2 = log10( rh(2:size(rh,1)-1,2:end-1) );
u2 = squeeze(ww(2,:,:)./ww(1,:,:)); u2 = u2(2:end-1,2:end-1);
v2 = squeeze(ww(3,:,:)./ww(1,:,:)); v2 = v2(2:end-1,2:end-1);
% subtract circular velocity
for jr=1:nr-2, v2(:,jr) = v2(:,jr)-aux.v0m(jr+1); end
% add symmetric part
rh2 = [rh2; rh2];
u2  = [u2 ; u2];
v2  = [v2 ; v2];
% also for patches
pat = reshape(pat,(nphi-1)*(nrr-1),4);
vrt = reshape(vrt,nphi*nrr,2);
% _____________________________________________________________________________
clf; 
if(0), 
  % as in Astron. Astrophys. 156 (1986), 354–380.
  cmg = colormap(gray); colormap(flipud(cmg)); cmapc = 1;
else
  % emphasize difference from axisymmetric case
  cmg = g2color; cmapc = 2;
end
% plot solution
mplot = [1 2 4];
for kplot=1:length(mplot),
  % select variable
  switch(kplot),
    case 1, h1 = rh2; ktitle = sprintf('log_{10} \\rho %s',tsc);
    case 2, h1 = u2;  ktitle = 'u';
    case 3, h1 = v2;  ktitle = 'v-v_0';
  end
  h1 = reshape(h1,numel(h1),1);
  % create plot
  subplot(2,2,mplot(kplot));
  patch('Faces',pat,'Vertices',vrt,'FaceVertexCData',h1,...
    'FaceColor','flat','EdgeColor','none','CDataMapping','scaled');
  axis equal tight; xlabel('x'); ylabel('y'); 
  if(cmapc == 2), ax = max(max(abs(h1))); caxis(ax*[-1 1]); end
  colorbar;
  title(ktitle);
end
%
% plot average density as a function of radius
subplot(4,2,5); semilogy(grd.r(2:end-1),rha(2:end-1));
xlabel('R'); ylabel('< \rho >_\phi');
%
% optional convergence history
if(nargin > 3 && length(conv) > 1),
  icmax = find(conv > 0);
  if(isempty(icmax)), icmax = length(conv); else icmax = icmax(end); end
  icmax = max(1,icmax);
  subplot(4,2,7); semilogy((0:length(conv)-1),conv/conv(1));
  xlabel('Iteration'); ylabel('Residual');
  axis([0 icmax 1.e-16 1.e2]);
end
% refresh plot
set(0,'CurrentFigure',myFigureHandle);
%EOF

