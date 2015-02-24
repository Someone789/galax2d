% makefigs.m
% prepare figures for a paper
% _____________________________________________________________________________
clear all; close all;
%
exepath = '../g2';
domirror = 1;
% select a file with results
% label = 'G07'; runpath = '.'; 
label   = input('Label used in g2run (without path): ','s');
runpath = input('Directory path (. for current): ','s');
if(length(runpath)>0 & runpath(end) ~= '/'), runpath = [runpath '/']; end
% _____________________________________________________________________________
fname = [runpath label 'w512.mat']; 
if(~exist(fname,'file')), 
  fname = [runpath label 'w256.mat']; 
  if(~exist(fname,'file')), fprintf(1,'%s not found\n',fname); return; end
end
% label for output *.eps files (hardCopy=0 for no files)
fnamepl = ['figbw' label '_']; hardCopy = 1;
% line style etc.
thick = 0.8;  % line thickness
%
ctype = 0; % colormap (0=gray,1=color)
% zoomed figures: size and spacing and clip
xmax = [-1 8]; dx = [0.04 0.01]; hclip = 0.05;
% _____________________________________________________________________________
% assuming that matlab files reside in exepath
eval(sprintf('addpath %s -end',exepath));
% _____________________________________________________________________________
load(fname); fprintf(1,'Loaded %s\n',fname);
if(xmax(1) <= 0), xmax(1) = gpar.Rmax; end
% old style
if(isfield(aux,'u2m')), aux.v0m = aux.u2m; aux.v0l = aux.u2l; end
% _____________________________________________________________________________
% overview as plotted in code
g2pltr(ww,grd,aux,conv,1);
if(0 & hardCopy), 
  figsave(fnamepl,1);
end
% _____________________________________________________________________________
% each plot separately
% log_10(rho),u,v-v0
nSmooth = 0; iscale = 1; [qq,rha] = g2wtoq(ww,grd,aux,iscale,nSmooth);
qq(1,:,:) = qq(1,:,:)/log(10); % change to log10
if(iscale)
  tsc = ' (scaled per radius)'; % text for title
else
  tsc = ''; % text for title
end
% patches
phi = (0:2*grd.nphi)*grd.dphi; rr = grd.rl(2:end-1);
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
% variables for plot; remove ghost points
qq = qq(:,2:end-1,2:end-1);
% add symmetric part
rh2 = squeeze(qq(1,:,:)); rh2 = [rh2; rh2];
u2  = squeeze(qq(2,:,:)); u2 = [u2; u2];
v2  = squeeze(qq(3,:,:)); v2 = [v2; v2];
% finalize patches
pat = reshape(pat,(nphi-1)*(nrr-1),4);
vrt = reshape(vrt,nphi*nrr,2);
% _____________________________________________________________________________
subplot(1,1,1); clf; 
switch(ctype),
  case 0, cmg = g2color(0); useColor = 0; % grey
  case 1, cmg = g2color(1); useColor = 1;  
end
% plot solution
titles = {sprintf('log_{10}\\rho%s',tsc),'u/c','(v-v_0)/c'};
ext = {'a','b','c'}; % for labelling figure files
kplot = 1;
if(0), % not used ...
  for k=1:3,
    % select variable
    switch(k),
      case 1, h1 = rh2;
      case 2, h1 = u2; 
      case 3, h1 = v2; 
    end
    h1 = reshape(h1,numel(h1),1); clf;
    % create plot
    patch('Faces',pat,'Vertices',vrt,'FaceVertexCData',h1,...
      'FaceColor','flat','EdgeColor','none','CDataMapping','scaled');
    axis equal tight; xlab='x'; ylab='y'; titl = titles{k}; 
    if(ctype == 1), h1x = max(max(abs(h1))); caxis(h1x*[-1 1]); end
    colorbar vert; dothick; drawnow; pause(0.1);
    if(hardCopy),
      fnamep = sprintf('%s%d%s',fnamepl,kplot,ext{k});
      figsave(fnamep,useColor);
    end
  end
end
% _____________________________________________________________________________
% plot average density as a function of radius
kplot = kplot+1; clf; 
% kplot=2;
if(0), % not used, see below, sequence
  semilogy(grd.r(2:end-1),rha(2:end-1),'k-');
  xlab ='R'; ylab = '<\rho>_\phi'; dothick; pause(0.1);
  if(hardCopy),
    fnamep = sprintf('%s%d',fnamepl,kplot); figsave(fnamep,0);
  end
end
% _____________________________________________________________________________
% convergence history
kplot = kplot+1;
icmax = find(conv > 0);
if(isempty(icmax)), icmax = length(conv); else icmax = icmax(end); end
icmax = max(1,icmax);
clf; semilogy((0:length(conv)-1),conv/conv(1),'k-');
xlab ='Iteration'; ylab = 'Residual'; 
axis([0 round(icmax*(1+1.e-3)) 1.e-16 1.e2]);
dothick; pause(0.1);
if(hardCopy), 
  fnamep = sprintf('%s%d',fnamepl,kplot); figsave(fnamep,0);
end
% _____________________________________________________________________________
% sequence of grids
kplot = kplot+1; iscale3 = 0;
ih = strfind(fname,'.mat'); % not failsafe
if(ih > 1), ih = ih-1; while(~isletter(fname(ih))), ih = ih-1; end; end
if(ih > 1),
  fname = fname(1:ih); nf = 0; fres = cell(10,2);
  for k=1:10,
    nk = 2^(2+k); fnamk = sprintf('%s%d.mat',fname,nk);
    if(exist(fnamk,'file')),
      nf = nf+1; fres{nf,1} = fnamk; fres{nf,2} = nk;%fprintf(1,'%s\n',fnamk);
    end
  end
  % read at most 3 files
  ksel = (max(1,nf-2):nf); nsel = length(ksel); 
  aruv = cell(nsel,4); qg = cell(nsel,1); ag = cell(nsel,1);
  ltype = {'k-.','k--','k-'};    %legn = cell(nsel,1);
  ltype = ltype(end+1-nsel:end); legn1 = {};
  %
  for k1=1:nsel, k = ksel(k1);
    fnamk = char(fres{k,1}); fprintf(1,'%s\n',fnamk);
    load(fnamk);
    qg{k1} = g2wtoq(ww,grd,aux,iscale3,0); % {log(rho),u,v-v0}, not log10
    % correct old style names
    if(isfield(aux,'u2m')), aux.v0m = aux.u2m; aux.v0l = aux.u2l; end
    ag{k1} = aux;  
    % 
    aruv{k1,1} = grd.r;
    for ic=1:3, 
      a = squeeze(qg{k1}(ic,2:end-1,:));
      if(ic == 1), a = exp(a); end; % rho
      aruv{k1,ic+1} = mean(a,1); % average over phi
    end
    legn1{k1} = sprintf('%d',fres{k,2});
  end
  % plot average of rho, u/c, (v-v_0)/c
  for mplot=1:3,
    switch(mplot),
      case 1, qcomp = 1; ylab = '<\rho>_\phi';      ftmp = 'rho';
      case 2, qcomp = 2; ylab = '<u/c>_\phi';       ftmp = 'u';
      case 3, qcomp = 3; ylab = '<(v-v_0)/c>_\phi'; ftmp = 'v1';
    end
    for k1=1:nsel,
      if(mplot == 1), 
        semilogy(aruv{k1,1},aruv{k1,mplot+1},ltype{k1}); 
      else
        plot(    aruv{k1,1},aruv{k1,mplot+1},ltype{k1});
      end
      if(k1==1), hold on; end
    end
    hold off; 
    xlabel('R'); legn = legn1; legloc = 'SouthEast'; 
    dothick; pause(0.1);
    if(hardCopy),
      fnamep = sprintf('%s%d%s',fnamepl,kplot,ftmp); figsave(fnamep,0);
    end
  end
end
% _____________________________________________________________________________
% zoomed images
ctype = 0; iscale2 = 0; titles{1} ='log_{10}\rho';
%
usensm = 1; % if 0, gaussian, not so good?
if(usensm == 0),
  sigmafac = 0.2; % no good
  sigmafac = 0.05;
else
  nsm = 128*[1,1,1]; % number of smoothing steps
  gs = []; 
end
clf; kplot = kplot+1; useColor = 0; g2color(0,2); 
for m=1:length(xmax),
  [qp,xp] = g2plotprepx(ww,grd,aux,iscale2,nSmooth,xmax(m),dx(m),domirror);
  dxp = mean(diff(xp)); 
  %
  if(usensm == 0),
    gsigma = sigmafac*xmax(m);
    fprintf(1,'  xmax = %g, dxp = %g, gsigma = %g\n',xmax(m),dxp,gsigma);
    % 1d conv in x and y
    hxp = exp(-xp.^2/(2*gsigma^2)); 
    ixp = find(hxp > 1.e-6); 
    xxp =  xp(ixp(1):ixp(end));
    hxp = hxp(ixp(1):ixp(end)); hxp = hxp/sum(hxp);
    if(1),
      gs = hxp'*hxp; % ?? zeros?
      if(1), imagesc(xxp,xxp,gs'); pause(0.1); end
      sgs = sum(sum(gs)); gs = gs/sgs;
      fprintf(1,'sum(gs)=%g -> 1\n',sgs);
      % or compare to nsm result? problem if sgs != 1?
    end
  else
    fprintf(1,'  xmax = %g, dxp = %g, nSmooth = %d,%d,%d\n',xmax(m),dxp,nsm);
  end
  %
  for k=1:3,
    % difference with smoothed version
    a = qp{k};
    if(usensm),
      for k1=1:nsm(k), a = smooth2(a); end
    else
      if(0), a = conv2(a,gs,'same'); % slow 
      else a = conv2(a,hxp,'same'); a = conv2(a,hxp','same'); end
    end
    a = qp{k}-a;
    imagesc(xp,xp,a'); colorbar vert; set(gca,'YDir','normal');
    amn = min(min(a)); amx = max(max(a));
    caxis(hclip*max(abs([amn,amx]))*[-1 1]); 
    % colormap(flipud(colormap()));
    axis equal tight; xlab ='x'; ylab = 'y';
    titl = sprintf('%s   (%d%% clip)',titles{k},round(100*hclip));
    legn = []; dothick; pause(0.1);
    % colormapeditor; 
    if(hardCopy),
      fnamep = sprintf('%s%dx%d%s%d',fnamepl,kplot,m,ext{k},usensm);
      figsave(fnamep,useColor);
    end
  end
end
% _____________________________________________________________________________
% performance
kplot = kplot+1;
p2 = perf(:,end)/60^2;
loglog(perf(:,1),p2,'ko-'); % axis tight;
axis([0.9*min(perf(:,1)) 1.1*max(perf(:,1)) 0.9*min(p2) 1.1*max(p2)]);
% xlab = 'n_R = n_\phi'; 
xlab = 'Number of points in R or \phi';
ylab = 'Cpu time (hours)';
dothick;
if(hardCopy),
  fnamep = sprintf('%s%d',fnamepl,kplot); figsave(fnamep,0);
end
%
if(size(perf,2)==3),
  kplot = kplot+1;
  loglog(perf(:,1),perf(:,2),'ko-'); axis tight;
  axis([0.9*min(perf(:,1)) 1.1*max(perf(:,1)) ...
        0.9*min(perf(:,2)) 1.1*max(perf(:,2))]);
  xlab = 'Number of points in R or \phi';
  ylab = 'Number of iterations';
  dothick;
  if(hardCopy),
    fnamep = sprintf('%s%d',fnamepl,kplot); figsave(fnamep,0);
  end
end
% _____________________________________________________________________________
%EOF




