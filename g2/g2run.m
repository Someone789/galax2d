% g2run.m for Matlab                            Copyright (c) 2013, W.A. Mulder
%
% - computes 2D stationary isothermal hydrodynamical galactic model
% - version with full Jacobian for 2nd order discretization and direct solver
% - please edit g2myparms.m to change the default parameters in g2setparm.m
% - if earlier results are found, an automatic restart is carried out
%
% _____________________________________________________________________________
clear all; close all; fclose all;
%
gpar = G2Parameters; % default parameters in class G2Parameters 
% _____________________________________________________________________________
% change parameters by using the file g2myparms.m, if any
% _____________________________________________________________________________
fchange = 'g2myparms.m';
if(exist(fchange,'file')),
  gpar0 = gpar; eval(strrep(fchange,'.m',''));
  ndif = g2pardif(gpar0,gpar); % print differences
  if(ndif > 0), 
    fprintf(1,'Changed default parameters using file %s\n',fchange);
  end
end
% _____________________________________________________________________________
% check parameters
gpar.ni = 2^ceil(log2(gpar.ni)); % force mesh size to be a power of 2
gpar.nf = 2^ceil(log2(gpar.nf));
gpar.resfactor2 = min(gpar.resfactor1,gpar.resfactor2);
gpar.order = max(1,min(2,gpar.order));
% _____________________________________________________________________________
% initialisation
% grid on finest mesh
grd = g2mesh(gpar.pp,gpar.kappa,gpar.Rmin,gpar.Rmax,gpar.nf,gpar.nf);
% gravitational potential
grv = g2pot2(gpar);
% _____________________________________________________________________________
grdf = grd;             % properties of final, finest grid
ni  = gpar.ni;          % initial grid has ni cells per coordinate
nf  = gpar.nf;          % final   grid has nf cells per coordinate
nlevel = 1+log2(nf/ni); % number of grid levels for successive refinement
kordermax = gpar.order; % 1st or 2nd order discretization
conv = [];              % store convergence history
perf = zeros(nlevel,3); % performance (number of points, step count, cpu time)
%
res = 0; res0 = 0; hasload = 0;
% _____________________________________________________________________________
% check if this is a restart
mlevel = 1; 
ctm = {'','t'};  % final or intermediate result
for klevel=nlevel:-1:1, % from coarse to find grids
  for k=1:2,
    fsave = sprintf('%sw%d%s.mat',gpar.label,ni*2^(klevel-1),ctm{k});
    if(exist(fsave,'file')),
      fprintf(1,'Found earlier result: %s\n',fsave);
      mlevel = klevel; % adjust starting level (grid size)
      break;
    end
  end
  if(mlevel > 1), break; end
end
% _____________________________________________________________________________
for klevel=mlevel:nlevel, % loop for successive refinement, coarse to fine
  tcpu0 = cputime();
  % subsample mesh
  if(klevel == mlevel),
    % modify gpar.nf
    gpar.nf = ni*2^(klevel-1); grd = g2meshcoarsen(grdf,gpar.nf);
    % initial solution and other auxiliary data
    [ww,aux] = g2init(gpar,grv,grd);
  end
  % order 1 or 2 or both?
  if(    gpar.nf <  gpar.norderswitch), korder1 = 1; korder2 = korder1;
  elseif(gpar.nf == gpar.norderswitch), korder1 = 1; korder2 = kordermax;
  else                                  korder1 = kordermax; korder2 = korder1;
  end
  %
  korder = korder1-1; kstep = 0;
  while(korder < korder2)
    korder = korder+1;
    if(gpar.print > -2), 
      fprintf(1,'--- level %d: n = %d, order = %d\n',klevel,grd.nrad,korder);
    end
    gpar.order = korder; % modify gpar.order
    if(gpar.print > 1), spparms('spumoni',2); end; % print sparse matrix info
    % filename for possible restart
    hasold = 0; nw = size(ww); 
    fsave = sprintf('%sw%d.mat',gpar.label,nw(2)-2);
    if(hasload == 0),
      if(exist(fsave,'file')),
        hasold = 1;
      else
        % filename for intermediate, not-yet-converged result (*t.mat)
        fsave = sprintf('%sw%dt.mat',gpar.label,nw(2)-2);
        if(exist(fsave,'file')), hasold = -1; end
      end
      if(hasold ~= 0), % restart from file
        load(fsave); 
        hasload = 1; fprintf(1,'Read solution from %s\n',fsave);
        % optional modification of parameters (not failsafe)
        fnewparms = [gpar.label 'newparms.m'];
        if(exist(fnewparms,'file')), 
          gpar0 = gpar; 
          fprintf(1,'Modifying parameters from restart file %s with %s\n',...
                  fsave,fnewparms);
          eval(strrep(fnewparms,'.m',''));
          ndif = g2pardif(gpar0,gpar); % print differences
          if(ndif > 0), hasold = -1; end; % treat as intermediate results
          fsave = sprintf('%sw%dt.mat',gpar.label,nw(2)-2);
          % this does not work with gpar.nf ... needs to be changed 
          % by g2myparms.m
        else
          fprintf(1,'%s not found, using original parameters\n',fnewparms);
        end
        % file has precedence
        if(gpar.order < korder),
          korder = gpar.order;
          fprintf(1,'Changed order to %d, given in file.\n',korder);
        end
      end
    end
    if(hasold == 1),
      continue; % skip iterations
    end
    % _________________________________________________________________________
    % iterations with Newton's method
    kbias = 0; krefr = 0;
    for k=1:gpar.nstep,
      if(korder == 2 && mod(k,gpar.nsave)==1),
        if(k > 1),
          fsave = sprintf('%sw%dt.mat',gpar.label,nw(2)-2);
          g2save; % save intermediate result
        end
        % recompute bias
        %   grd = g2bias(ww,grd,aux); fprintf(1,'    Bias recomputed\n');
        % simpler: just take a small value
        grd.bphi=1.e-10+0*grd.bphi; grd.brad=1.e-10+0*grd.brad; kbias = k;
      end
      %
      % rhs and matrix
      rh = g2rhs(ww,1,aux,gpar,grd);
      % local factor for timestep and residual norm
      [edt,res] = g2tim1(ww,rh.rhs,grd);
      % convergence history
      conv = [conv;res];
      % plot latest
      g2pltr(ww,grd,aux,conv,1); pause(0.1);
      %
      if(k == 1), % store first value
        if(hasold > -1), res0 = res;
        else res0 = max(conv); 
        end
      else
        if(gpar.print > -1), fprintf(1,'%5d: %g\n',k,res/res0); end
        if(res > 1.e4*res0),
          fprintf(1,'Divergence?\n'); break;
        elseif( (korder < 2 && res < gpar.resfactor1*res0) | ...
                              (res < gpar.resfactor2*res0))
          break;         
        end;
      end
      % solve linear system and update solution
      [ww,bfac,relch] = g2sol(ww,rh,gpar.idtfactor*edt,gpar.relchange,...
                              gpar.print);
      spparms('spumoni',0); % suppress further output
      % clean up solution
      ww = g2wfix(ww,aux,gpar,grd);
      kstep = kstep+1;
    end
    % _________________________________________________________________________
  end; % korder
  % ___________________________________________________________________________
  if(hasold ~= 1),
    dcpu = cputime-tcpu0;
    if((size(perf,1) < klevel) || (perf(klevel,2) == 0)),
      perf(klevel,:) = [gpar.nf kstep dcpu];
    elseif(perf(klevel,1) == gpar.nf) 
      perf(klevel,3) =  perf(klevel,2)+dcpu;
    else fprintf(1,'perf?\n');
    end
    if(gpar.print > -2), 
      fprintf(1,'    %d steps, %g s\n',kstep,perf(klevel,3)); 
    end
    % divergence (probably)
    if(res > 1.e4*res0), break; end
    % save
    if((gpar.nf > gpar.norderswitch && gpar.order > 1) || (gpar.nf > 128)),
      fsave = sprintf('%sw%d.mat',gpar.label,nw(2)-2);
      g2save;
    end
  end
  hasold = 0;
  % ___________________________________________________________________________
  % refine solution, doubling the number of grid cells per coordinate
  if(klevel < nlevel),
    grdc = grd; auxc = aux;
    gpar.nf = 2*gpar.nf; % double number of grid cells per coordinate
    % map final finest grid to coarser grid
    grd = g2meshcoarsen(grdf,gpar.nf);
    % auxiliary data on finer grid
    [w0,aux] = g2init(gpar,grv,grd); clear w0;
    % interpolate to finer mesh
    ww = g2refine(ww,gpar, grdc,auxc, grd,aux);
    % plot result after interpolation
    g2pltr(ww,grd,aux,[],1);
  end
end
% _____________________________________________________________________________
%EOF

