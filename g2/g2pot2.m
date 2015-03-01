function grv = g2pot2(gpar)
% function grv = g2pot2(gpar)
% creates struct with values for gravitational potential
%
% input: struct gpar with pp,axs,axi,ii,om
%    pp       power for density as function of radius (e.g., -1.8)
%    axs      ratio for short        axis/long axis (axs < axi < 1)
%    axi      ratio for intermediate axis/long axis
%    ii       power for cut-off function (e.g., 10)
%    om       rotation  (e.g., 0.1)
%    cutoff   1: bar cut-off at co-rotation (CR),
%             2: cut-off at the Outer Lindblad Resonance (OLR)
%    c        sound speed

% rho = rho_0 m^p, m^2=(x/a)^2+(y/b)^2+(z/c)^2, a = 1,
% note: axs = c/a, axi = b/a
% rho(r,theta,phi)=rho_0 r^p sum_{n,m} a_{nm}P_n^m(\mu)cos(m\phi),
% \mu=cos(theta), a_{nm}=2 as_{nm} for m !=0 and a_{n0}=as_{n0},
% n and m even, n >= m

% _____________________________________________________________________________
% compute a_{nm}, 0<=m<=n, both even, n <= nmax
nmax = 10; % maximum number of terms in series expansion, spherical harmonics
dotst = 0; % set to 1 for testing
% _____________________________________________________________________________
% longest axis a in x, shortest c in z
% c < b < a, a=1, ep1 = 1-(c/a)^2, ep2 = 1-(c/b)^2, ep = 1-ep2/ep1
% m^2=(r/c)^2*(1-h); h = ep1*(1-mu^2)*(1-ep*sin(phi)^2));
% axs = c/a, axi = b/a
ep1 = 1-gpar.axs^2; ep2 = 1-(gpar.axs/gpar.axi)^2; ep = 1-ep2/ep1;
pp = gpar.pp;
%
bm = zeros(nmax/2+1,1); bm(1) = 1;
for m=0:2:nmax-2,
  jm = 1+m/2; bm(jm+1) = bm(jm)*(-ep/4)/(1+m/2);
end
%
bn = zeros(nmax/2+1,1); bn(1) = gpar.axs^(-gpar.pp);
for n=0:2:nmax-2,
  jn = 1+n/2; bn(jn+1) = bn(jn)*(pp-n)*(ep1/2)/((2*n+3)*(2*n+1));
end
% _____________________________________________________________________________
% fill bnm and anm
bnm = zeros(nmax/2+1,nmax/2+1); anm = bnm;
for n=0:2:nmax, jn = 1+n/2;
  for m=0:2:n,  jm = 1+m/2;
    % b_{nm}
    fr = factorial(n-m)/factorial((n-m)/2); if(m > 0), fr = 2*fr; end
    bnm(jn,jm) = bn(jn)*bm(jm)*fr;
    % a_{nm}=b_{nm} sum_{k=0}^{\infty} t_k S_k
    % tsm = sum td, td = t_k S_k
    k = 0; tk = 1; Sk = sksum(k,n,m,ep); tsk = tk*Sk; tsm = tsk;
    % fprintf(1,'%d %d,%d: %g %g  %g\n',n,m,k, Sk,tk,tsm);
    while( abs(tsk) > 1.e-16*abs(tsm) )
      k = k+1;
      Sk = sksum(k,n,m,ep);
      ft = ((pp-n)/2+1-k)/(2*(k+n)+1) * (k+(n+m)/2)/k * (-2*ep1);
      tk = tk*ft;
      % test Sk and tk
      if(dotst && k < 20), % test
        jmax = k+(n-m)/2;
        jj = (0:jmax); sh1 = factorial(jmax)*factorial(m/2);
        sh2 = factorial(2*jj+m).*(-ep/4).^jj;
        sh3 = factorial(jmax-jj).*factorial(jj+m/2).* ...
              factorial(jj+m).*factorial(jj);
        Sk2 = sh1*sum(sh2./sh3);
        %
        tk2 = fracf(pp/2,k+n/2)/fracf(pp/2,n/2)*...
          factorial(k+(n+m)/2)/factorial((n+m)/2)*...
          factorial(k+n+1)/factorial(n+1)*...
          factorial(2*n+2)/factorial(2*k+2*n+2)*((-4*ep1)^k)/ ...
          factorial(k);
        srat = Sk/Sk2; trat = tk/tk2;
        if(abs(srat-1) > 1.e-12 || abs(trat-1) > 1.e-12),
          fprintf(1,'%d %d,%d: %g,%g|%g\t%g,%g|%g\n',n,m,k,Sk,Sk2, ...
            srat,tk,tk2,trat);
        end
      end
      tsk = tk*Sk; tsm = tsm+tsk;
      % fprintf(1,'%d %d,%d: %g %g  %g\n',n,m,k, Sk,tk,tsm);
    end
    anm(jn,jm) = bnm(jn,jm)*tsm;
  end
end
% _____________________________________________________________________________
% typical values:
% -1.8                 pp
%  0.5                 axs
%  0.8                 axi
%  0.1                 om
%  10                  ii
%  0.5570251e+00       a00
% -0.3735773e+00       a20
%  0.3760795e-01       a22
%  0.1394757e+00       a40
% _____________________________________________________________________________
pp = gpar.pp; ii = gpar.ii; om = gpar.om;
% c_{nm} = a_{nm}/( (p+2-n)*(p+3+n) ) if no cut-off
% P_n^m(\mu) for \mu=cos(theta)=0:
%    (-1)^(n/2)*(n+m)! / (2^n ((n+m)/2)! ((n-m)/2)! )
nmax = 4; % maximum value of n (even)
% axisymmetric contribution to potential, m=0
% NOTE: c00 is not the same as cnm in paper, but
% the sum of contributions cn0
c00 = 0;
for nh=0:nmax/2,
  c00 = c00+anm(nh+1,1)*Pnm0(2*nh,0)/((pp+2-2*nh)*(pp+3+2*nh) );
end
%
fa0 = sqrt(c00*(pp+2)); rcr = (om/fa0)^(2/pp); a22 = anm(2,2); % (1+n/2,1+m/2)
b = zeros(1,5); b(2) = 3.*a22/(pp*(pp+5)); % c22=a22/(p*(p+5)), P22(0)=3
% rotating bar with cut-off
if( (gpar.cutoff > 0) & (om ~= 0) ),
  switch( gpar.cutoff ),
    case 1, rcp = (om/fa0)^2; % cut-off at co-rotation
    case 2, fa1 = sqrt(pp+4); % cut-off at the Outer Lindblad Resonance
      rcp = (om/(fa0*(1+0.5*fa1)))^2;
  end
  rcut = rcp^(1/pp); rci = rcut^ii; h0 = -(6./5.)*a22*ii^2;
  b(1) =  h0*rcp/(pp*(pp+ii)*(pp+2*ii));
  b(3) = -6.*a22/((pp+  ii)*(pp+5+  ii)*rci  );
  b(4) =  3.*a22/((pp+2*ii)*(pp+5+2*ii)*rci^2);
  b(5) =  h0*rcp*rcut^5 /((pp+5)*(pp+5+ii)*(pp+5+2*ii));
else
  rcut = 1.e32;
end
% rescale potential by 1/c^2 (c is unit of velocity)
bc = b/gpar.c^2;
% output
grv = struct('a22',a22,'c00',c00,'fa0',fa0,'rcr',rcr,'rcut',rcut,'bc',bc);
% _____________________________________________________________________________
% initial values, circular velocity
% analytical solution at inner boundary -- see g2init.m
% potential and derivatives
% _____________________________________________________________________________
function a = Pnm0(n,m)
% function a = Pnm0(n,m)
% associated Legendre function at zero argument for even 0<=m<=n:
%   (-1)^k * (n+m)! /(2^n * k! * ((n-m)/2)!) with k=(n+m)/2
k = (n+m)/2;
a = (-1)^k * factorial(n+m) /( 2^n * factorial(k) * factorial((n-m)/2) );
if(1), % test
  b = legendre(n,0); b = b(m+1);
  if(abs(a/b-1) > 1.e-12),
    fprintf(1,'Problem in Pnm0(%d,%d)=%g != %g\n',n,m,b,a);
  end
end
% _____________________________________________________________________________
function Sk = sksum(k,n,m,ep)
% Sk = sum_{j=0}^{k+(n-m)/2} s_j for even 0<=m<=n
% eq.(AII.7), Astron. Astrophys. 134, 158--170 (1984)
Sk = 1; sj = 1; jmax = k+(n-m)/2;
for j=1:jmax,
  sj = sj*(-ep/2)*(2*j+m-1)/(j+m)*(jmax+1-j)/j; Sk = Sk+sj;
end
% _____________________________________________________________________________
function a = fracf(q,m)
% function a = fracf(q,m)
%  a = q*(q-1)*...*(q-(m-1)), m > 0;
%  a = 1 for m=0
a = 1;
if(m > 0),
  a = q;
  if(m > 1),
    b = q; for j=2:m, b = b-1; a = a*b; end;
  end
end
%EOF

