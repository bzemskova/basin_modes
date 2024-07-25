function params = gendata_params()

% simulation parameters
params.f = 1e-4; % coriolis parameter
params.deltaT = 100; % timestep
params.n_periods = 16;
params.g = 9.81;
params.Lz = 4000;
params.z_ratio = 4000/4000;
params.Lz_tot = 4000;

% wavenumber and forcing frequency
params.om = linspace(1.05,3,60)*params.f;
params.k = params.om/sqrt(params.g*params.Lz);

params.Lxbasin = 400e3; %basin width
params.alpha = 1; %basin aspect ratio
params.inlet_frac = 0.2; %inlet fraction

% filename to save basin coordinates
params.filename = 'basin_400km_400km.mat';