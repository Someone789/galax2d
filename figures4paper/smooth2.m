function b = smooth2(a)
% SMOOTH2	smooth 2d array a by convolution with [1/4 1/2 1/4] in
%               each coordinate
% b = smooth2(a)
nn = size(a);
ii = (1:nn(1)); im = ii-1; im(1)=1; ip = ii+1; ip(nn(1)) = nn(1);
jj = (1:nn(2)); jm = jj-1; jm(1)=1; jp = jj+1; jp(nn(2)) = nn(2);
b = 0.25*a+0.125*(a(im,jj)+a(ip,jj)+a(ii,jm)+a(ii,jp))+...
     0.0625*(a(im,jm)+a(im,jp)+a(ip,jm)+a(ip,jp));
%%EOF