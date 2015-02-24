function wf = g2refine(wc,gpar,grdc,auxc,grdf,auxf)
% wf = g2refine(wc,grdc,auxc,grdf,auxf)
% interpolate coarse-grid solution "wc" on grid "grdc" with auxiliaries "auxc"
% to "wf" on grid "grdf" with auxiliary variables "auxf"

wc = g2wfix(wc,auxc,gpar,grdc); % boundary conditions
nwc = size(wc); 
wc(2,:,:) = wc(2,:,:)./wc(1,:,:); % u
wc(3,:,:) = wc(3,:,:)./wc(1,:,:); % v
for jrc=1:nwc(2),
  wc(1,:,jrc) = wc(1,:,jrc)./grdc.r(jrc);  % rho
  wc(3,:,jrc) = wc(3,:,jrc)-auxc.v0m(jrc); % v-v0
end
% nwc(1) contains number of components and does not change
% nwc(2) and nwc(3) have 2 ghost cells
nwf = nwc; nwf(2:3) = 2*(nwf(2:3)-2)+2;
wf = zeros(nwf);
% interpolate rho,u,v-v0
for jrc=2:nwc(2)-1,
  jrf = 2*(jrc-1);
  for jpc=2:nwc(2)-1,
    jpf = 2*(jpc-1);
    wlc=0.5625*wc(:,jpc  ,jrc  );
    wml=0.1875*wc(:,jpc-1,jrc  );
    wmr=0.1875*wc(:,jpc+1,jrc  );
    wmb=0.1875*wc(:,jpc  ,jrc-1);
    wmt=0.1875*wc(:,jpc  ,jrc+1);
    wf(:,jpf  ,jrf  )=wlc+wml+wmb+0.0625*wc(:,jpc-1,jrc-1);
    wf(:,jpf+1,jrf  )=wlc+wmr+wmb+0.0625*wc(:,jpc+1,jrc-1);
    wf(:,jpf  ,jrf+1)=wlc+wml+wmt+0.0625*wc(:,jpc-1,jrc+1);
    wf(:,jpf+1,jrf+1)=wlc+wmr+wmt+0.0625*wc(:,jpc+1,jrc+1);
  end
end
% back to conserved variables
for jrf=1:nwf(2),
  wf(1,:,jrf) = wf(1,:,jrf).*grdf.r(jrf);
  wf(3,:,jrf) = wf(3,:,jrf)+auxf.v0m(jrf);
end
wf(2,:,:) = wf(2,:,:).*wf(1,:,:);
wf(3,:,:) = wf(3,:,:).*wf(1,:,:);
%
wf = g2wfix(wf,auxf,gpar,grdf); % boundary conditions
%EOF

