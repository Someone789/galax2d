function [A,b,nb] = g2mksa(ww,rh,idt)
% function [A,b,nb] = g2mksa(ww,rh,idt)
% build sparse matrix A from rh

if(isempty(rh.A00)), A = []; return; end;
if(nargin > 1 && ~isempty(idt)), hasidt = 1; else hasidt = 0; end
nw = size(rh.rhs); nn = nw-[0 2 2]; na = prod(nn); 
n1 = nn(1); n2 = nn(2); n3 = nn(3); n1q = n1^2;
% 1 main block diagonals and 2*4 off-diagonals, block has size n1q;
% naz can be sharpened!
naz = 9*n1^2*n2*n3; ia = zeros(1,naz); ja = ia; va = ia;
% if dorw: scale by 1/<R*rho>_R
dorw = 1; rw = 1+zeros(1,n3);
if(dorw),
  for j3=1:n3, rw(j3) = 1/mean(ww(1,:,j3+1)); end
end
%
rr3 = reshape(repmat((1:n1)',1,n1) ,1,n1q); %123123123
rc3 = reshape(repmat((1:n1)',1,n1)',1,n1q); %111222333
% cell (j2,j3) has matrix A**(:,:,j2+1,j3+1)
jj = 0; % count element nr of A
ddj = 0;
for j3=1:n3, % radial
  for j2=1:n2, % phi is periodic
    if(hasidt), ddj = idt(j2+1,j3+1)*eye(3); end; % 1/dt damping
    k0 = n1*((j2-1)+n2*(j3-1)); % row offset
    % each row has 5*n1^2 entries except near boundary
    % (j2,j3-2)B0m|(j2,j3-1)A0m|(j2-2,j3)Bm0|(j2-1,j3)Am0|
    %   (j2,j3)A00|
    % (j2+1,j3)Ap0|(j2+2,j3)Bp0|(j2,j3+1)A0p|(j2,j3+2)B0p
    %
    if(j3 > 1), % (j2,j3-1)A0m R
      l0 = k0-n1*n2; % column offset
      kn = jj+(1:n1q); % add n1^2 elements
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.A0m(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end  
    if(j3 > 2), % (j2,j3-2)B0m R
      l0 = k0-2*n1*n2; % column offset
      kn = jj+(1:n1q); % add n1^2 elements
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.B0m(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
    if(1), % (j2-1,j3)Am0 periodic phi
      l0 = k0-n1; % column offset
      if(j2 == 1), l0 = l0+n1*n2; end % periodic
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.Am0(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
    if(1), % (j2-2,j3)Bm0 periodic phi
      l0 = k0-2*n1; % column offset
      if(j2 <= 2), l0 = l0+n1*n2; end % periodic
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.Bm0(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
    if(1) % (j2,j3)A00 -- central
      l0 = k0; % column offset        
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = (rh.A00(:,:,j2+1,j3+1)-ddj)*rw(j3);
      jj = jj+n1q;
    end
    if(1) % (j2+1,j3)Ap0 periodic phi
      l0 = k0+n1; % column offset        
      if(j2 == n2), l0 = l0-n1*n2; end % periodic
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.Ap0(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
    if(1), % (j2+2,j3)Bp0 periodic phi
      l0 = k0+2*n1; % column offset        
      if(j2 >= n2-1), l0 = l0-n1*n2; end % periodic
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.Bp0(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
    if(j3 < n3), % (j2,j3+1)A0p R
      l0 = k0+n1*n2; % column offset
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.A0p(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end  
    if(j3<n3-1), % (j2,j3+2)B0p R
      l0 = k0+2*n1*n2; % column offset
      kn = jj+(1:n1q); 
      ia(kn) = k0+rr3; % rows
      ja(kn) = l0+rc3; % cols
      va(kn) = rh.B0p(:,:,j2+1,j3+1)*rw(j3);
      jj = jj+n1q;
    end
  end;
end;
%
ip = find(ia > 0); 
if(~isempty(ip)), ia = ia(ip); ja = ja(ip); va = va(ip); end
ip = find(ja > 0);
if(~isempty(ip)), ia = ia(ip); ja = ja(ip); va = va(ip); end
%
A = sparse(ia,ja,va,na,na); % contains A
%
b = rh.rhs(:,2:end-1,2:end-1);
for j3=1:n3, b(:,:,j3) =  b(:,:,j3)*rw(j3); end
nb = size(b);
b = reshape(b,numel(b),1); % make 1-D vector
%EOF
