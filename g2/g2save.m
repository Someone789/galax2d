% g2save.m
% script to save results to the file "fsave"
% used in g2run.m; sloppy coding ...
if(~exist('perf','var')), perf = []; end
save(fsave,'gpar','grd','grv','aux','ww','conv','res0','perf');
fprintf(1,'Results saved in %s\n',fsave);
%
