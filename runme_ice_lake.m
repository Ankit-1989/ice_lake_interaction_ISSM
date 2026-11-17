%% ========================================================================
%  runme_ice_lake.m
%
%  ISSM experiments for:
%    Pramanik, A., Greenwood, S. L., Morlighem, M., Barnett, J., Holmes, F. A.,
%    Gyllencreutz, R. and Regnell, C.: Modelling the effects of proglacial lake
%    filling and drainage on lake-terminating glacier dynamics.
%
%  Domain : Qooqqup Sermia catchment (Narsarsuaq, south Greenland)
%  Physics: SSA, Budd sliding law (p = q = 1), constant SMB (MAR mean),
%           no frontal ablation, prescribed lake level.
%
%  ------------------------------------------------------------------------
%  HOW THE LAKE IS REPRESENTED
%  The lake level is prescribed through ISSM's sea level field
%  (md.initialization.sealevel), but only on the nodes inside the lake
%  contour (Exp/lake_front_contour_points.exp). Sea level stays at 0 m
%  everywhere else. The lake floor is flattened to -200 m (Step 3), and
%  lake levels are given as WATER DEPTH above the lake floor, as in the paper:
%
%        sealevel (lake nodes) = depth + lake_floor     (e.g. 100 m -> -100 m)
%
%  Basal water pressure then follows Pw = max(0, rho_w g (sl - b)).
%
%  ------------------------------------------------------------------------
%  STEPS
%   Base model
%     1  Mesh               Anisotropic mesh refined on fast flow, fjord and lake
%     2  Parameterization   Geometry, masks, velocities, rheology (ice_lake.par)
%     3  Geometry_update    SSA, ice mask, flat lake floor (-200 m), BCs
%     4  Inversion_Budd     Basal friction coefficient (Budd) inversion
%   Spin-ups (moving front, von Mises calving, 20 yr)
%     5  SpinUp_Lake100     100 m lake
%     6  SpinUp_NoLake      no lake (0 m)
%     7  SpinUp_Lake200     200 m lake (floating terminus, start of Exp 3)
%   Experiments (fixed front)
%     8  Exp1_Lake100       Exp 1   : 100 m lake, 100 yr
%     9  Exp1_NoLake        Exp 1   : land-terminating, 100 yr
%    10  Exp2a_Rise         Exp 2a  : 100 -> 200 m in 50 yr
%    11  Exp2a_Fall         Exp 2a  : 100 ->   0 m in 50 yr
%    12  Exp2b_Rise         Exp 2b  : 0 -> 100/125/150/175/200 m in 50 yr
%    13  Exp2b_Fall         Exp 2b  : 100/125/150/175/200 -> 0 m in 50 yr
%    14  Exp3a_Drainage     Exp 3a  : 100 m drop in 1/2/3/4 yr, from 200 m and 100 m
%    15  Exp3b_GLOF         Exp 3b  : 10/20/30/40/50 m drop in 1 day, from 200 m
%
%  USAGE
%    steps = 1:4;   runme_ice_lake     % build and invert the base model
%    steps = 5:7;   runme_ice_lake     % spin-ups
%    steps = 8:15;  runme_ice_lake     % experiments
%  Each step loads the output of the step it depends on from ./Outputs.
%
%  REQUIREMENTS
%    ISSM (https://github.com/ISSMteam/ISSM), MATLAB R2016b or newer.
%    Data   : Model_Data/vx.nc, Model_Data/vy.nc  (MEaSUREs velocity mosaic)
%             Model_Data/BedMachineGreenland-v5.nc
%             Model_Data/lake_mask.mat
%             Model_Data/MARv3.11.3-historical-combined.nc  (monthly SMB)
%             Exp/NSQ_mouginot_modified_basin_w_LdFjLk.exp (domain outline)
%             Exp/lake_front_contour_points.exp            (lake area)
%    Helper functions (not part of ISSM's core):
%             interpBedmachineGreenland, interpMAR_monthly, interpSeaRISE_new
% =========================================================================

%% ======================== User settings =================================
if ~exist('steps','var'), steps = 1:4; end

addpath('./Scripts');   % folder with the helper functions listed above

% Domain and mesh
glacier           = 'Narsarsuaq';
exp_file          = './Exp/NSQ_mouginot_modified_basin_w_LdFjLk.exp';
parameterize_file = './ice_lake.par';
bedmachine        = './Model_Data/BedMachineGreenland-v5.nc';
hmin      = 300;     % minimum mesh size (m)
hmax      = 1000;    % maximum mesh size (m)
fjordmesh = 100;     % target mesh size in fjord and lake (m)

% Lake setup
cfg.lake_contour = './Exp/lake_front_contour_points.exp';
cfg.lake_floor   = -200;   % flattened lake floor elevation (m)
cfg.ocean_sl     = 0;      % sea level outside the lake (m)
cfg.rho_water    = 1000;   % fresh lake water density (kg m^-3)

% Spin-up front / calving (von Mises, spin-ups only)
cfg.sigma_grounded        = 2e6;    % Pa
cfg.sigma_floating        = 150e3;  % Pa (ISSM default)
cfg.min_calving_thickness = 100;    % m

% SMB: MAR monthly records averaged to a constant field (1960-2000 mean).
% Adjust the record indices to the time axis of your MAR file.
cfg.smb_file   = './Model_Data/MARv3.11.3-historical-combined.nc';
cfg.smb_months = 1:40*12;
cfg.bedmachine = bedmachine;

% Time stepping
cfg.dt      = 1/12;       % model time step (yr), monthly
cfg.outfreq = 1/cfg.dt;   % save output once a year
nyrs_spinup = 20;         % spin-up length (yr)
nyrs_exp    = 100;        % experiment length (yr)

% Computing
cfg.np  = 6;      % number of MPI processes
doplot  = true;   % quick-look figures for the base steps

org = organizer('repository','Outputs','prefix',glacier,'steps',steps);
cfg.outdir = 'Outputs';
cfg.prefix = glacier;

%% ======================== BASE MODEL ====================================

%% Step 1: Mesh
if perform(org,'Mesh')

    md = bamg(model,'domain',exp_file,'hmax',2000);

    % Observed velocities, clipped to the domain
    x  = double(ncread('Model_Data/vx.nc','x'));
    y  = double(ncread('Model_Data/vx.nc','y'));
    vx = double(ncread('Model_Data/vx.nc','Band1'));
    vy = double(ncread('Model_Data/vy.nc','Band1'));

    posx  = find(x<=max(md.mesh.x) & x>=min(md.mesh.x));
    posy  = find(y> min(md.mesh.y) & y<=max(md.mesh.y));
    x_dom = x(posx);
    y_dom = y(posy);
    vx_dom = vx(posx,posy)';
    vy_dom = vy(posx,posy)';

    md.inversion.vx_obs  = InterpFromGridToMesh(x_dom,y_dom,vx_dom,md.mesh.x,md.mesh.y,0);
    md.inversion.vy_obs  = InterpFromGridToMesh(x_dom,y_dom,vy_dom,md.mesh.x,md.mesh.y,0);
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);

    % Refine where ice is fast (> 500 m/yr), in the fjord and in the lake.
    % BedMachine mask: 0 = ocean, 1 = ice-free land (lakes are mapped as land).
    hmaxVertices = NaN*ones(md.mesh.numberofvertices,1);
    hmaxVertices(md.inversion.vel_obs>500) = hmin;
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest',bedmachine);
    hmaxVertices(M==0 | M==1) = fjordmesh;

    md = bamg(md,'hmin',hmin,'hmax',hmax,'field',md.inversion.vel_obs,...
        'gradation',1.1,'err',25,'anisomax',1.,'hmaxVertices',hmaxVertices);

    % Re-interpolate velocities on the final mesh
    md.inversion.vx_obs  = InterpFromGridToMesh(x_dom,y_dom,vx_dom,md.mesh.x,md.mesh.y,0);
    md.inversion.vy_obs  = InterpFromGridToMesh(x_dom,y_dom,vy_dom,md.mesh.x,md.mesh.y,0);
    md.inversion.vel_obs = sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);

    if doplot
        plotmodel(md,'data',md.inversion.vel_obs,'title','Observed velocity (m/yr)',...
            'data','mesh','title','Mesh');
    end
    savemodel(org,md);
end

%% Step 2: Parameterization
if perform(org,'Parameterization')

    md = loadmodel(org,'Mesh');
    md = parameterize(md,parameterize_file);

    if doplot
        plotmodel(md,'data',md.mask.ocean_levelset,'title','Ocean mask',...
            'data',md.geometry.thickness,'title','Thickness (m)',...
            'data',md.geometry.bed,'title','Bed (m)',...
            'data','BC','caxis#1',[-1 1]);
    end
    savemodel(org,md);
end

%% Step 3: Geometry update and forward stress balance
if perform(org,'Geometry_update')

    md = loadmodel(org,'Parameterization');
    md = setflowequation(md,'SSA','all');
    md.levelset.stabilization = 1;

    % Ice mask from BedMachine (mask >= 2: ice), as a signed distance
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest',bedmachine);
    md.mask.ice_levelset(M<2) = +1;
    md.mask.ice_levelset(md.mask.ice_levelset==0) = -1;   % no zeros (advection)
    md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    % Flat lake floor: remove the influence of lake bathymetry
    pos = find(md.mesh.y<-3.178e6 & md.mesh.x<0.15e4 & md.geometry.bed<-150);
    md.geometry.bed(pos)  = cfg.lake_floor;
    md.geometry.thickness = max(md.geometry.surface-md.geometry.base,0);

    md.masstransport.spcthickness = NaN(md.mesh.numberofvertices,1);

    % Boundary conditions: no flow across ice-covered domain boundary,
    % observed velocities on ice-free land
    pos = find(md.mesh.vertexonboundary & M>=2);
    md.stressbalance.spcvx(pos) = 0;
    md.stressbalance.spcvy(pos) = 0;
    pos = find(M==1 | ~isnan(md.stressbalance.spcvx));
    md.stressbalance.spcvx(pos) = md.inversion.vx_obs(pos);
    md.stressbalance.spcvy(pos) = md.inversion.vy_obs(pos);

    md.inversion.iscontrol = 0;
    md = solve(md,'sb');
    savemodel(org,md);
end

%% Step 4: Inversion of the Budd friction coefficient
%  tau_b = C^2 N u_b  (p = q = 1), N = rho_i g H + rho_w g b
if perform(org,'Inversion_Budd')

    md = loadmodel(org,'Geometry_update');

    md.inversion = m1qn3inversion();
    md.inversion.vx_obs  = md.initialization.vx;
    md.inversion.vy_obs  = md.initialization.vy;
    md.inversion.vz_obs  = md.initialization.vz;
    md.inversion.vel_obs = md.initialization.vel;

    md.inversion.iscontrol = 1;
    md.inversion.maxsteps  = 70;
    md.inversion.maxiter   = 70;

    % Cost: absolute misfit (101) + log misfit (103) + Tikhonov regularisation (501)
    md.inversion.cost_functions = [101 103 501];
    md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3);
    md.inversion.cost_functions_coefficients(:,1) = 2000;
    md.inversion.cost_functions_coefficients(:,2) = 1;
    md.inversion.cost_functions_coefficients(:,3) = 6e-7;

    md.inversion.control_parameters = {'FrictionCoefficient'};
    md.inversion.min_parameters = 10*ones(md.mesh.numberofvertices,1);
    md.inversion.max_parameters = 300*ones(md.mesh.numberofvertices,1);

    % Rheology: ice at -5 C
    md.materials.rheology_B(:) = cuffey(273.15-5);

    % Ice mask: everything north of y = -3.183e6 m is ice; shift the
    % ice-free mask by one element to avoid a cliff at the transition
    md.mask.ice_levelset(md.mesh.y>-3.183e6) = -1;
    pos = find(max(md.mask.ice_levelset(md.mesh.elements),[],2)>0);
    md.mask.ice_levelset(md.mesh.elements(pos,:)) = 1;
    md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);

    md.stressbalance.restol = 0.01;
    md.stressbalance.reltol = 0.1;
    md.stressbalance.abstol = NaN;

    % Budd friction
    md.friction.p = ones(md.mesh.numberofelements,1);
    md.friction.q = ones(md.mesh.numberofelements,1);
    md.friction.coupling = 2;

    % No damage, no thermal model
    md.damage.D         = zeros(md.mesh.numberofvertices,1);
    md.damage.spcdamage = NaN*ones(md.mesh.numberofvertices,1);
    md.thermal.spctemperature = NaN*ones(md.mesh.numberofvertices,1);

    md.toolkits = toolkits();
    md.verbose  = verbose('solution',false,'module',false,'convergence',false);
    md.cluster  = generic('name',oshostname,'np',cfg.np);
    md = solve(md,'sb');

    % Keep inverted friction; replace values on floating/ice-free/deep-bed
    % areas with a bed-dependent value (Akesson et al., 2018)
    md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;
    pos = find(md.mask.ocean_levelset<0 | md.mask.ice_levelset>1 | ...
        (md.friction.coefficient>150 & md.mask.ocean_levelset<1000 & md.geometry.bed<0));
    md.friction.coefficient(pos) = 120*(min(max(0,md.geometry.bed(pos)+800),max(md.geometry.bed))/max(md.geometry.bed));

    % Geometry consistency: base not below bed, base = bed where grounded
    pos = find(md.geometry.bed>md.geometry.base);
    md.geometry.base(pos) = md.geometry.bed(pos);
    pos = find(md.mask.ocean_levelset>0);
    md.geometry.base(pos) = md.geometry.bed(pos);
    md.geometry.surface   = md.geometry.base+md.geometry.thickness;

    % Boundary conditions for the transient runs
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest',bedmachine);
    pos = find(md.mesh.vertexonboundary & M>=1);
    md.stressbalance.spcvx(pos) = 0;
    md.stressbalance.spcvy(pos) = 0;
    pos = find(M==1 & md.mesh.vertexonboundary==0);
    md.stressbalance.spcvx(pos) = NaN;
    md.stressbalance.spcvy(pos) = NaN;

    if doplot
        plotmodel(md,'data',md.inversion.vel_obs,'title','Observed velocity (m/yr)',...
            'data',md.results.StressbalanceSolution.Vel,'title','Modelled velocity (m/yr)',...
            'data',md.friction.coefficient,'title','Friction coefficient',...
            'data',md.results.StressbalanceSolution.Vel-md.inversion.vel_obs,'title','Modelled - observed (m/yr)',...
            'caxis#1-2',[1.5 round(max(md.inversion.vel_obs),-2)],'log#1-2',10,...
            'colorbar#all','on','ncols',4);
    end
    savemodel(org,md);
end

%% ======================== SPIN-UPS ======================================
%  Moving front with von Mises calving, constant lake level, 20 years.

%% Step 5: Spin-up with a 100 m lake
if perform(org,'SpinUp_Lake100')
    md = setup_spinup(loadmodel(org,'Inversion_Budd'),cfg);
    md = set_lake_level(md,100,[],cfg);
    md = set_timestepping(md,cfg.dt,nyrs_spinup,cfg.outfreq);
    md = run_transient(md,cfg);
    savemodel(org,md);
end

%% Step 6: Spin-up without a lake (land-terminating)
if perform(org,'SpinUp_NoLake')
    md = setup_spinup(loadmodel(org,'Inversion_Budd'),cfg);
    md = set_lake_level(md,0,[],cfg);
    md = set_timestepping(md,cfg.dt,nyrs_spinup,cfg.outfreq);
    md = run_transient(md,cfg);
    savemodel(org,md);
end

%% Step 7: Spin-up with a 200 m lake (terminus floating)
if perform(org,'SpinUp_Lake200')
    md = setup_spinup(loadmodel(org,'Inversion_Budd'),cfg);
    md = set_lake_level(md,200,[],cfg);
    md = set_timestepping(md,cfg.dt,nyrs_spinup,cfg.outfreq);
    md = run_transient(md,cfg);
    savemodel(org,md);
end

%% ======================== EXPERIMENTS ===================================
%  All experiments use a fixed ice front and restart from a spin-up.
%  Lake level is given as [depths] at [times] (yr from the run start) and is
%  interpolated linearly in between and held constant after the last time.

%% Step 8: Exp 1 - lake-terminating, fixed 100 m lake
if perform(org,'Exp1_Lake100')
    md = restart_at(loadmodel(org,'SpinUp_Lake100'));
    md = set_lake_level(md,100,[],cfg);
    md = run_experiment(md,cfg,nyrs_exp);
    savemodel(org,md);
end

%% Step 9: Exp 1 - land-terminating (no lake)
if perform(org,'Exp1_NoLake')
    md = restart_at(loadmodel(org,'SpinUp_NoLake'));
    md = set_lake_level(md,0,[],cfg);
    md = run_experiment(md,cfg,nyrs_exp);
    savemodel(org,md);
end

%% Step 10: Exp 2a - lake level rise, 100 -> 200 m over 50 yr
if perform(org,'Exp2a_Rise')
    md = restart_at(loadmodel(org,'SpinUp_Lake100'));
    md = set_lake_level(md,[100 200],[0 50],cfg);
    md = run_experiment(md,cfg,nyrs_exp);
    savemodel(org,md);
end

%% Step 11: Exp 2a - lake level fall, 100 -> 0 m over 50 yr
if perform(org,'Exp2a_Fall')
    md = restart_at(loadmodel(org,'SpinUp_Lake100'));
    md = set_lake_level(md,[100 0],[0 50],cfg);
    md = run_experiment(md,cfg,nyrs_exp);
    savemodel(org,md);
end

%% Step 12: Exp 2b - lake level rise from 0 m to D over 50 yr
exp2b_depths = [200 175 150 125 100];
if perform(org,'Exp2b_Rise')
    md_spinup = loadmodel(org,'SpinUp_NoLake');
    for D = exp2b_depths
        md = restart_at(md_spinup);
        md = set_lake_level(md,[0 D],[0 50],cfg);
        md = run_experiment(md,cfg,nyrs_exp);
        save_run(md,cfg,sprintf('Exp2b_Rise_0to%d',D));
    end
end

%% Step 13: Exp 2b - lake level fall from D to 0 m over 50 yr
%  Starts from year 70 of the matching Exp 2b rising run (lake stable at D),
%  holds D for 1 yr, then drains linearly to 0 m over 50 yr.
t_restart_2b = 70;
if perform(org,'Exp2b_Fall')
    for D = exp2b_depths
        md = restart_at(load_run(cfg,sprintf('Exp2b_Rise_0to%d',D)),t_restart_2b);
        md = set_lake_level(md,[D D 0],[0 1 51],cfg);
        md = run_experiment(md,cfg,nyrs_exp);
        save_run(md,cfg,sprintf('Exp2b_Fall_%dto0',D));
    end
end

%% Step 14: Exp 3a - 100 m drainage in 1, 2, 3 and 4 years
%  From a floating terminus (200 m lake) and a grounded terminus (100 m lake).
%  Lake is held for t_hold_3a years, then drops by 100 m over N years.
drain_years = [1 2 3 4];
t_hold_3a   = 1;     % yr before drainage starts
nyrs_3a     = 10;    % run length (yr)
if perform(org,'Exp3a_Drainage')
    start = struct('depth',{200,100},'spinup',{'SpinUp_Lake200','SpinUp_Lake100'});
    for s = start
        md_spinup = loadmodel(org,s.spinup);
        for N = drain_years
            md = restart_at(md_spinup);
            md = set_lake_level(md,[s.depth s.depth s.depth-100],[0 t_hold_3a t_hold_3a+N],cfg);
            md = run_experiment(md,cfg,nyrs_3a);
            save_run(md,cfg,sprintf('Exp3a_from%d_drop100_in%dyr',s.depth,N));
        end
    end
end

%% Step 15: Exp 3b - GLOF-like drainage of 10-50 m in one day
%  From a floating terminus (200 m lake). Hourly time step; the lake is held
%  for t_hold_3b days, drops linearly over 24 hours, and the run continues
%  to 90 days.
glof_drops   = [10 20 30 40 50];   % m
dt_3b        = 1/365/24;           % hourly time step (yr)
outfreq_3b   = 24;                 % output every 24 steps (daily); 1 = hourly
t_hold_3b    = 1;                  % days before drainage starts
ndays_3b     = 90;                 % run length (days)
if perform(org,'Exp3b_GLOF')
    md_spinup = loadmodel(org,'SpinUp_Lake200');
    cfg3b = cfg;
    cfg3b.dt = dt_3b;
    cfg3b.outfreq = outfreq_3b;
    for dd = glof_drops
        md = restart_at(md_spinup);
        md = set_lake_level(md,[200 200 200-dd],[0 t_hold_3b t_hold_3b+1]/365,cfg);
        md = run_experiment(md,cfg3b,ndays_3b/365);
        save_run(md,cfg,sprintf('Exp3b_GLOF_drop%dm_1day',dd));
    end
end

%% ======================== Local functions ===============================

function md = setup_spinup(md,cfg)
% Common transient settings for the spin-ups: constant SMB, no basal melt,
% grounding-line migration, moving front with von Mises calving.

    md.inversion.iscontrol = 0;
    md.materials.rho_water = cfg.rho_water;

    % No damage, no thermal model, no hydrology
    md.damage.D         = zeros(md.mesh.numberofvertices,1);
    md.damage.spcdamage = NaN*ones(md.mesh.numberofvertices,1);
    md.transient.isthermal    = 0;
    md.transient.ishydrology  = 0;
    md.thermal.spctemperature = NaN*ones(md.mesh.numberofvertices,1);

    % SMB: mean of the monthly MAR fields (mm w.e./month -> m ice eq./yr)
    disp('Reading and interpolating MAR SMB');
    smb = zeros(md.mesh.numberofvertices,numel(cfg.smb_months));
    for k = 1:numel(cfg.smb_months)
        smb(:,k) = interpMAR_monthly(md.mesh.x,md.mesh.y,'SMB',cfg.smb_months(k),cfg.smb_file);
    end
    smb(smb==-9999) = 0;
    md.smb.mass_balance = mean(smb,2)/1000*12*(md.materials.rho_freshwater/md.materials.rho_ice);

    % No basal melt
    md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);

    % Geometry consistency
    pos = find(md.geometry.bed>md.geometry.base);
    md.geometry.base(pos) = md.geometry.bed(pos);
    md.geometry.surface   = md.geometry.base+md.geometry.thickness;

    % Grounding line
    md.transient.isgroundingline = 1;
    md.groundingline.migration              = 'SubelementMigration';
    md.groundingline.melt_interpolation     = 'NoMeltOnPartiallyFloating';
    md.groundingline.friction_interpolation = 'SubelementFriction1';

    % Moving front: reset ice levelset from BedMachine, fixed on the boundary
    md.transient.ismovingfront = 1;
    M = interpBedmachineGreenland(md.mesh.x,md.mesh.y,'mask','nearest',cfg.bedmachine);
    md.mask.ice_levelset(M<2) = +1;
    md.mask.ice_levelset(md.mask.ice_levelset==0) = -1;
    md.mask.ice_levelset = reinitializelevelset(md,md.mask.ice_levelset);
    md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
    bnd = find(md.mesh.vertexonboundary);
    md.levelset.spclevelset(bnd) = md.mask.ice_levelset(bnd);
    md.levelset.kill_icebergs = 1;

    % von Mises calving, no frontal melt
    md.calving = calvingvonmises();
    md.calving.stress_threshold_groundedice = cfg.sigma_grounded;
    md.calving.stress_threshold_floatingice = cfg.sigma_floating;
    md.calving.min_thickness = cfg.min_calving_thickness;
    md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);

    md.transient.requested_outputs = {'TotalSmb','SmbMassBalance',...
        'IceVolume','IceVolumeAboveFloatation','IceMass',...
        'IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled',...
        'GroundedArea','FloatingArea','MaskOceanLevelset','MaskIceLevelset',...
        'TotalFloatingBmb','TotalGroundedBmb','BasalforcingsFloatingiceMeltingRate',...
        'GroundinglineMassFlux','IcefrontMassFluxLevelset',...
        'TotalCalvingFluxLevelset','CalvingMeltingrate','TotalCalvingMeltingFluxLevelset',...
        'Calvingratex','Calvingratey','CalvingCalvingrate','SigmaVM',...
        'DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy',...
        'DeviatoricStressxz','DeviatoricStressyz','DeviatoricStresszz'};
end

function md = set_lake_level(md,depths,times,cfg)
% Prescribe the lake level (water depth above the lake floor, m) on the lake
% nodes. depths: scalar (constant level) or vector with matching times (yr).

    lake = find(ContourToNodes(md.mesh.x,md.mesh.y,cfg.lake_contour,2));
    sl = cfg.ocean_sl*ones(md.mesh.numberofvertices,numel(depths));
    sl(lake,:) = repmat(depths(:)'+cfg.lake_floor,numel(lake),1);

    if numel(depths)==1
        md.initialization.sealevel = sl;
    else
        md.initialization.sealevel = [sl; times(:)'];
    end
end

function md = set_timestepping(md,dt,nyrs,outfreq)
% Fixed time step, run starts at t = 0.

    md.timestepping = timestepping();
    md.timestepping.time_step  = dt;
    md.timestepping.start_time = 0;
    md.timestepping.final_time = nyrs;
    md.settings.output_frequency = outfreq;
end

function md = run_experiment(md,cfg,nyrs)
% Fixed-front transient run used by all experiments.

    md = set_timestepping(md,cfg.dt,nyrs,cfg.outfreq);
    md.transient.ismovingfront = 0;
    md.transient.requested_outputs = {'TotalSmb','SmbMassBalance',...
        'IceVolume','IceVolumeAboveFloatation','IceMass',...
        'IceVolumeAboveFloatationScaled','GroundedAreaScaled','FloatingAreaScaled',...
        'GroundedArea','FloatingArea','MaskOceanLevelset','MaskIceLevelset',...
        'TotalFloatingBmb','BasalforcingsFloatingiceMeltingRate',...
        'GroundinglineMassFlux','IcefrontMassFluxLevelset','SigmaVM',...
        'DeviatoricStressxx','DeviatoricStressyy','DeviatoricStressxy',...
        'DeviatoricStressxz','DeviatoricStressyz','DeviatoricStresszz'};
    md = run_transient(md,cfg);
end

function md = run_transient(md,cfg)
% Solve the transient and warn if the run stopped early.

    md.toolkits = toolkits();
    md.toolkits.DefaultAnalysis = bcgslbjacobioptions();
    md.cluster = generic('name',oshostname,'np',cfg.np);
    md.verbose = verbose('solution',true,'module',false,'convergence',false);

    md = solve(md,'Transient');

    t_end = md.results.TransientSolution(end).time;
    if abs(t_end-md.timestepping.final_time) > md.timestepping.time_step
        warning('Run stopped at t = %g yr (final time %g yr).',t_end,md.timestepping.final_time);
    end
end

function md = restart_at(md,t)
% Restart from a previous transient run at time t (default: last output).

    sol = md.results.TransientSolution;
    if nargin<2
        idx = numel(sol);
    else
        [~,idx] = min(abs([sol.time]-t));
    end
    md = transientrestart(md,idx);
    md.results = struct();
end

function save_run(md,cfg,tag)
% Save one run of a multi-run experiment as Outputs/<prefix>_<tag>.mat
    save(fullfile(cfg.outdir,[cfg.prefix '_' tag '.mat']),'md','-v7.3');
end

function md = load_run(cfg,tag)
% Load a run saved with save_run
    S  = load(fullfile(cfg.outdir,[cfg.prefix '_' tag '.mat']),'md');
    md = S.md;
end
