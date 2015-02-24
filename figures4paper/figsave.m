function figsave(filename,useColor)
% function figsave(filename)
% save current figure to filename.eps and filename.fig
if(length(filename) < 1), return; end
if(nargin > 1 && useColor), 
  ptype = 'epsc2'; 
else 
  ptype = 'eps2';
end
cm = sprintf('print -d%s %s.eps',ptype,filename); eval(cm);
cm = sprintf('print -dtiff %s.tif',filename); eval(cm);
savefig([filename '.fig']);
%
fprintf(1,'Plot saved as %s.eps and %s.tif and %s.fig\n',...
  filename,filename,filename);
%EOF
