function [wn,bfac,rch] = g2sol(ww,rh,idt,relchangemax,doprint)
% function[wn,bfac,rch] = g2sol(ww,rh,idt,relchangemax)
% solve  A (wn-ww) = -rhs(ww) for wn
% - fill sparse matrix: A = -idt+d(rhs)/dw, negative definite
% - include idt if not empty
% - solve and update solution with a maximum relative change
%   of relchangemax
% - extra printout of doprint > 0

% input: relchangemax,doprint; idt may be empty
% target change of rchmax < 1
if(nargin > 3 && relchangemax > 0), rchmax = relchangemax; 
else rchmax = 0.6; end
if(nargin > 4), doprint1 = doprint; else doprint1 = -1; end
% _____________________________________________________________________________
% fill matrix and right-hand side: A = -idt+d(rhs)/dw, b = rhs
[A,b,nb] = g2mksa(ww,rh,idt);
n1 = nb(1); n2 = nb(2); n3 = nb(3); 
% b is interior of rh.rhs
if(isempty(b)),
  b = rh.rhs(:,2:end-1,2:end-1);
  nb = size(b); n1 = nb(1); n2 = nb(2); n3 = nb(3);
  b = reshape(b,numel(b),1); % make 1D vector
end
% _____________________________________________________________________________
% solve
b = A\b; bminmax = [min(b),max(b)];
b = reshape(b,n1,n2,n3);
% measure relative change, assuming c=1, so rho*c=ww(1,:,:)
rch = [0 0 0]; jp = (2:n2+1); jr = (2:n3+1);
ch1 = b(1,jp-1,jr-1)./( ww(1,jp,jr) );
ch2 = b(2,jp-1,jr-1)./( ww(1,jp,jr) + abs(ww(2,jp,jr)) );
ch3 = b(3,jp-1,jr-1)./( ww(1,jp,jr) + abs(ww(3,jp,jr)) );
cha = [max(max(abs(ch1))) max(max(abs(ch2))) max(max(abs(ch3)))];
for k=1:3, rch(k) = max(rch(k),cha(k)); end
%
rch2 = max(rch); 
if(rch2 > rchmax),
  bfac = rchmax/rch2;
else
  bfac = 1;
end
% update: wn = ww-A^{-1}rhs 
wn = zeros(size(ww)); 
wn(:,2:end-1,2:end-1) = ww(:,2:end-1,2:end-1)-bfac*b;
%
if(doprint1 > 0),
  fprintf(1,'change in w: [%g,%g], relative: %g %g %g',bminmax,rch);
  if(bfac ~= 1), fprintf(1,'reduced by %g',bfac); end
  fprintf(1,'\n');
end
%EOF
