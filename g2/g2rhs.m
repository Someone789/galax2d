function rh = g2rhs(ww,doA, aux,gpar,grd)
% function rh = g2rhs(ww,doA, aux,gpar,grd)
%
% given conserved states ww, compute struct rh with fields rhs
% if doA=1, the jacobian A = d rhs/d w is computed as well
% A consists of
%          B0p
%          A0p
%  Bm0 Am0 A00 Ap0 Bp0
%          A0m
%          B0m
%
% aux,gpar,grd contain global variables

ww = g2wfix(ww,aux,gpar,grd); % boundary conditions
% initialize
nw = size(ww); rhs = zeros(3,nw(2),nw(3));
if(doA), A00 = zeros(3,3,nw(2),nw(3)); else A00 = [ ]; end
Am0 = A00; Ap0 = A00; A0m = A00; A0p = A00;
Bm0 = A00; Bp0 = A00; B0m = A00; B0p = A00;
%
omc = gpar.om/gpar.c; % omega/c
oc2 = 2*omc;
% _____________________________________________________________________________
% w = R rho {1,u,v}, f = R rho {u,u^2+c^2,uv}, g = rho {v,uv,v^2+c^2}
% dw/dt=s-df/dR-dg/d\phi
% s = rho {0, -R dV/dR + c^2 + (v+om R)^2,-dV/d\phi-u(v+2 om R)}
% _____________________________________________________________________________
% source term s times local volume dR d\phi
for jr=2:nw(3)-1, % source term, radial coordinate
  for jp=2:nw(2)-1,
    %rho = ww(1,jp,jr)/grd.r(jr);     % w1/Rm
    uu  = ww(2,jp,jr)/ww(1,jp,jr); % u
    vv  = ww(3,jp,jr)/ww(1,jp,jr); % v
    vr  = vv/grd.r(jr); % v/R
    h21 = aux.b21(jp,jr)-vv*vr; h23 =  (oc2+2*vr);
    h31 = aux.b31(jp,jr)+uu*vr; h32 = -(oc2+  vr); h33 = -uu/grd.r(jr);
    As = [0 0 0;h21 0 h23;h31 h32 h33]*grd.dvol(jr);
    rhs(:,jp,jr) = As*ww(1:3,jp,jr);
    if(doA), A00(:,:,jp,jr) = As; end
  end
end
% _____________________________________________________________________________
% flux g in phi: g = rho {v,uv,v^2+c^2}
for jr=2:nw(3)-1, % R
  dp = g2gradphi(ww(:,:,jr),grd,jr,gpar.order,doA);
  dw = dp.druv;
  %
  for jp=2:nw(2)-1, % phi
    rho = dp.rho(jp); uu = dp.u(jp); vv = dp.v(jp);      
    jpl = jp-1; if(jpl <       2), jpl = jpl+(nw(2)-2); end; % periodic
    jpr = jp+1; if(jpr > nw(2)-1), jpr = jpr-(nw(2)-2); end; % periodic
    % left  side of cell
    vb = vv-dw(3,jp);
    if(vb < 1)
      %fprintf(1,'phi,L: %d %d\n',jp,jpl);
      rhb = rho-dw(1,jp);
      ub  = uu -dw(2,jp);
      % flux
      [g,Ag] = g2fphi(-1,grd.dr(jr),rhb,ub,vb);
      rhs(:,jp ,jr)=rhs(:,jp ,jr)+g;
      rhs(:,jpl,jr)=rhs(:,jpl,jr)-g; % jp to the right of jpl
      if(doA)
        AL =         -Ag*dp.BL(:,:,jp) *dp.A1(:,:,jpl); 
        AM =  Ag*(eye(3)-dp.BM(:,:,jp))*dp.A1(:,:,jp ); 
        AR =         -Ag*dp.BR(:,:,jp) *dp.A1(:,:,jpr); 
        Am0(:,:,jp ,jr)=Am0(:,:,jp ,jr)+AL;
        A00(:,:,jp ,jr)=A00(:,:,jp ,jr)+AM;
        Ap0(:,:,jp ,jr)=Ap0(:,:,jp ,jr)+AR;
        A00(:,:,jpl,jr)=A00(:,:,jpl,jr)-AL;
        Ap0(:,:,jpl,jr)=Ap0(:,:,jpl,jr)-AM;
        Bp0(:,:,jpl,jr)=Bp0(:,:,jpl,jr)-AR;
      end
    end
    % right side of cell
    vb = vv+dw(3,jp);
    if(vb > -1)
      rhb = rho+dw(1,jp);
      ub  = uu +dw(2,jp);
      % flux
      [g,Ag] = g2fphi(+1,grd.dr(jr),rhb,ub,vb);
      rhs(:,jp ,jr)=rhs(:,jp ,jr)-g;
      rhs(:,jpr,jr)=rhs(:,jpr,jr)+g; % jp to the left of jpr
      if(doA)
        AL =          Ag*dp.BL(:,:,jp) *dp.A1(:,:,jpl); 
        AM =  Ag*(eye(3)+dp.BM(:,:,jp))*dp.A1(:,:,jp ); 
        AR =          Ag*dp.BR(:,:,jp) *dp.A1(:,:,jpr); 
        Am0(:,:,jp ,jr)=Am0(:,:,jp ,jr)-AL;
        A00(:,:,jp ,jr)=A00(:,:,jp ,jr)-AM;
        Ap0(:,:,jp ,jr)=Ap0(:,:,jp ,jr)-AR;
        Bm0(:,:,jpr,jr)=Bm0(:,:,jpr,jr)+AL;
        Am0(:,:,jpr,jr)=Am0(:,:,jpr,jr)+AM;
        A00(:,:,jpr,jr)=A00(:,:,jpr,jr)+AR;
      end
    end
  end; % jp (phi)
end;
% _____________________________________________________________________________
% flux f in R
% note that f also has interpolation if 1st order
% _____________________________________________________________________________
for jp=2:nw(2)-1, % phi
  rhr = squeeze(ww(1,jp,:))./grd.r';     % w1/Rm
  uur = squeeze(ww(2,jp,:)./ww(1,jp,:));
  vvr = squeeze(ww(3,jp,:)./ww(1,jp,:));
  v1r = vvr-aux.v0m'; % relative to axisymm. case
  % half gradients per cell
  dp = g2gradr(squeeze(ww(:,jp,:)),grd,aux,gpar.order,doA); 
  dw = dp.druv;
  % R, left  side of cell, jr=nw(3) values have to be reset
  for jr=2:nw(3), % note: 1 extra cell for Rmax boundary
    uu = uur(jr); ub = uu-dw(2,jr);
    if(ub <  1)
      jrl = jr-1; jrr = jr+1;
      rho = rhr(jr); v1 = v1r(jr);
      rhb = rho-dw(1,jr);
      vb  = v1 +aux.v0l(jr)-dw(3,jr);
      % radial flux
      [f,Af] = g2frad(-1,grd.rl(jr),grd.dphi,rhb,ub,vb);
      rhs(:,jp,jr )=rhs(:,jp,jr )+f;
      rhs(:,jp,jrl)=rhs(:,jp,jrl)-f; % jr to the right of jrl
      if(doA && jr < nw(3)), % at left side of cell
        AL =         -Af*dp.BL(:,:,jr) *dp.A1(:,:,jrl); 
        AM =  Af*(eye(3)-dp.BM(:,:,jr))*dp.A1(:,:,jr ); 
        AR =         -Af*dp.BR(:,:,jr) *dp.A1(:,:,jrr); 
        A0m(:,:,jp ,jr)=A0m(:,:,jp ,jr)+AL;
        A00(:,:,jp ,jr)=A00(:,:,jp ,jr)+AM;
        A0p(:,:,jp ,jr)=A0p(:,:,jp ,jr)+AR;
        A00(:,:,jp,jrl)=A00(:,:,jp,jrl)-AL;
        A0p(:,:,jp,jrl)=A0p(:,:,jp,jrl)-AM;
        B0p(:,:,jp,jrl)=B0p(:,:,jp,jrl)-AR;
      end
    end
  end
  % R, right side of cell, jr=1 values have to be reset  
  for jr=1:nw(3)-1, % note: 1 extra for Rmin boundary
    uu = uur(jr); ub = uu+dw(2,jr);
    if(ub > -1)
      jrl = jr-1; jrr = jr+1;
      rho = rhr(jr); v1 = v1r(jr);
      rhb = rho+dw(1,jr); 
      vb  = v1 +aux.v0l(jr+1)+dw(3,jr);
      % radial flux
      [f,Af] = g2frad(+1,grd.rl(jr+1),grd.dphi,rhb,ub,vb);
      rhs(:,jp,jr )=rhs(:,jp,jr )-f;
      rhs(:,jp,jrr)=rhs(:,jp,jrr)+f;
      if(doA)
        if(jrl > 0), AL = Af*dp.BL(:,:,jr) *dp.A1(:,:,jrl); 
        else         AL = zeros(3,3); 
        end
        AM =  Af*(eye(3)+dp.BM(:,:,jr))*dp.A1(:,:,jr ); 
        AR =          Af*dp.BR(:,:,jr) *dp.A1(:,:,jrr); 
        A0m(:,:,jp ,jr)=A0m(:,:,jp ,jr)-AL;
        A00(:,:,jp ,jr)=A00(:,:,jp ,jr)-AM;
        A0p(:,:,jp ,jr)=A0p(:,:,jp ,jr)-AR;
        B0m(:,:,jp,jrr)=B0m(:,:,jp,jrr)+AL; % used if jrr>3
        A0m(:,:,jp,jrr)=A0m(:,:,jp,jrr)+AM; % used if jrr>2
        A00(:,:,jp,jrr)=A00(:,:,jp,jrr)+AR;
      end
    end
  end; % jr (R)
end; % jp (phi)
% 
A0m(:,:,:,      2) = 0; % lower R boundary
B0m(:,:,:,      2) = 0; 
B0m(:,:,:,      3) = 0;
A0p(:,:,:,nw(2)-1) = 0; % upper R boundary
B0p(:,:,:,nw(2)-1) = 0; 
B0p(:,:,:,nw(2)-2) = 0;
% _____________________________________________________________________________
rh = struct('rhs',rhs,'A00',A00,...
            'Am0',Am0,'Bm0',Bm0,'Ap0',Ap0,'Bp0',Bp0,...
            'A0m',A0m,'B0m',B0m,'A0p',A0p,'B0p',B0p);
%EOF
