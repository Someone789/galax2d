% g2test.m                                      Copyright (c) 2013, W.A. Mulder
% runs one example
% please edit g2myparms.m to choose your own parameters
% _____________________________________________________________________________
% path to source code
addpath '../g2' -end;
%
if(0),
  % maximum number of processors
  maxNumCompThreads(4); % Matlab complains but does not offer an alternative?
  fprintf(1,'_____________________________________________________________________________\n');
end
% start
g2run;
%EOF

