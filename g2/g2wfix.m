function ww = g2wfix(ww,aux,gpar,grd)
% function ww = g2wfix(ww,aux,gpar,grd)
% - checks for negative densities in ww
% - imposes boundary values of ww(1:3,:,1) at rmin and of ww(1:3,:,end) at rmax
%   using info from aux and grd
% - imposes periodicity in the phi direction

% check for negative density
nw = size(ww); ineg = 0;
for jr=2:nw(3)-1,
  for jp=2:nw(2)-1,
    if(ww(1,jp,jr) <= 0),
      fprintf(1,'---- ww(1,%d,%d) = %g;\n',jp,jr,ww(1,jp,jr));
      ineg = ineg+1; ww(1,jp,jr) = 1.e-32;
    end
  end
end
if(ineg > 0),
  fprintf(1,'---- fixed zero or negative density in %d cells\n',ineg);
end
%
% set outer boundary
ww(1,:,end) = gpar.rhoouter*grd.r(end);
ww(2,:,end) = 0;
ww(3,:,end) = ww(1,:,end)*aux.v0l(end);
%
% set inner boundary
rhc = gpar.rhoinner;
ww(1,:,1) = aux.cen(1,:).*rhc;
ww(2,:,1) = aux.cen(2,:).*ww(1,:,1);
ww(3,:,1) = aux.cen(3,:).*ww(1,:,1);
% adjust periodicity (again)
ww = g2wperiod(ww);
%EOF

