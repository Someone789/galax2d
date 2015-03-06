% G2Parameters.m                                Copyright (c) 2013, W.A. Mulder
% class holding user and code parameters, defaults are set here
% does not work with Octave, use g2setpar() instead
% _____________________________________________________________________________
% Filenames:
%    label          for naming output files
%
% Physical parameters:
%    c              sound speed, in model units
%                   velocities u and v will have c as unit
%    rhoinit        initial density in interior
%    rhoinner       average density at inner boundary with radius Rmin [>>1]

%    rhoouter       average density at outer boundary with radius Rmax [1]
%    Rmin           minimum radius of domain
%    Rmax           maximum radius of domain
%
% Gravitational potential:
%   background density responsible for gravitational field:
%   rho_g = rho_{g,0} m^p, with m^2=(x/a)^2+(y/b)^2+(z/c)^2, a = 1,
%   note: axs = c/a, axi = b/a
%   approximated by
%      \rho_g(r,theta,phi)=\rho_{g,0} r^p sum_{n,m} a_{nm}P_n^m(\mu)cos(m\phi),
%      \mu=cos(\theta), n and m even, n >= m
%   ----------------------------------------------------------------------------
%    pp             power for gravitational density as a function
%                   of radius [-1.8]
%    axs            ratio of short to long axis (axs < axi < 1),
%                   the shortest axis points out of plane in the z-direction
%    axi            ratio of intermediate to long axis
%                   the longest axis points in the x-direction
%    ii             power for cut-off function
%    om             frame rotation rate, equals 2 pi/T, T is period
%    cutoff         0: no cutuff
%                   1: bar cut-off at co-rotation,
%                   2: cut-off at the Outer Lindblad Resonance
%    
% Numerical parameters:
%    ni             initial mesh has ni^2 cells for domain
%                      [Rmin,Rmax]x[0,pi]      
%    nf             final mesh   has nf^2 cells 
%                      ni and nf should be powers of 2
%                   the code performs successive grid refinement
%    kappa          R grid is equidistant in xi = r^(1+kappa*pp/2) 
%    order          spatial discretization of order 1 or 2
% 
% Parameters related to the iterative scheme:
%    idtfactor      scaling of 1/dt factor 
%                   Too small a value leads to smaller solution changes
%                   and may slow down the convergence.
%                   Too a large value may push to solution too far
%                   away from the stationary case.
%    relchange      max. relative change in solution per step
%    nstep          max. number of iterations per successive refinement level
%    nsave          iteration interval between saving intermediate results
%    norderswitch   switch from 1st to 2nd when n = norderswitch, n
%                   being the number of cells per coordinate
%    resfactor1     desired residual reduction for 1st order scheme
%    resfactor2     desired residual reduction for 2nd order scheme
%    print          extra printout for values >= 0
% _____________________________________________________________________________
classdef G2Parameters
  properties
    % for naming output files
    label     = 'G01';
    % physics
    c         = 0.035; % sound speed
    rhoinit   =     1; % initial density in interior
    rhoinner  =  1.e2; % average density at inner boundary
    rhoouter  =     1; % average density at outer boundary
    Rmin      =  0.25; % minimum radius of domain
    Rmax      =    30; % maximum radius of domain
    % gravitational potential
    pp     = -1.8; % power for gravitational density as a function
                   % of radius [-1.8]
    axs    =  0.5; % ratio of short to long axis (axs < axi < 1),
    axi    =  0.8; % ratio of intermediate to long axis 
    om     =  0.1; % frame rotation rate (2 pi / period)
    ii     =   10; % power for cut-off function for bar
    cutoff =    1; % 0:no cutuff
                   % 1: bar cut-off at co-rotation,
                   % 2: cut-off at the Outer Lindblad Resonance
    % numerical parameters
    ni        =     8; % initial mesh has ni^2 cells for domain
    nf        =   256; % final   mesh has nf^2 cells
    kappa     =     1; % R grid is equidistant in xi = r^(1+kappa*pp/2)
    order     =     2; % order of spatial discretization
    % iterative scheme            
    idtfactor =     2; % scaling of 1/dt factor
    relchange =   0.9; % max. relative change in solution per step
    nstep     =  4000; % max. number of iterations per level
    nsave     =    50; % interval between saving intermediate results
    norderswitch =    64;  % switch order 1 to 2 at norderswitch^2 cells
    resfactor1   = 1.e-8;  % residual reduction for 1st order scheme
    resfactor2   = 1.e-12; % residual reduction for 2nd order scheme
    %
    print = -1;            % extra printout for values >= 0
  end
end
%EOF


