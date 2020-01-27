
%% Description
%  This script sets all parameters that define the model and inversion
%  for the true inversion with the real data

%% Inversion parms
inv = struct(    'verbose',false                 ,... % option to spit out more information+plots
                 'niter',1000                   ,... % Number of iterations
                 'burnin',7000                   ,... % don't record results before burnin iterations
                 'cooloff',4000                  ,... % # of iterations over which temperature declines as erf
                 'tempmax',5                     ,... % maximum multiple of all standard deviations
                 'saveperN',30                   ,... % save only every saveperN iterations       
                 'bestNmod2keep',-5000           ,... % keep only the best N models in each chain, defined here
                 'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
                 'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'maxnkchain',350                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'nchains',10                    ,... % number of chains to start in parallel
                 'Nsavestate',2500               ,... % Niter per which the state of the parallel inversion is saved in .mat file
                 'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
                 'BWclust',1                     ,... % option to use only one c x             
                 'datatypes',{{'RF_Ps','RF_Sp','SW_Ray_phV','SW_Lov_phV'}});  
                                % any of {{'SW_x_y' with x='Ray/Lov' and y='phV/grV'; 
                                %          'BW_x_y' with x='Sp/Ps' and y=' /lo/fl';}}
                                %          'RF_x_y' with x='Sp/Ps' and y=' /CCP';}}
                                %          'HKstack_x' with x='P'

%% Model parms
modl = struct([]);

modl(1).nstas = 1;
modl.maxz = 300;                                      % maximum depth in model from ref ellipsoid, km
modl.maxkz = 250;                                     % maximum depth of deepest non-basal knot, km
modl.dz = 2;                                        % depth spacing of model, km

modl.sed = struct(...
    ... thickness of the sediments
                     'hmax',0.0                  ,... %5 max sed layer thickness, km
                     'hmin',0.0                  ,... %0 min sed layer thickness, km
                     'hstd',0.0                  ,... % std of sed layer thickness for perturbation, km
    ... shear velocity of the sediments
                     'vsmax',3.3                 ,... % max sed velocity, km/s
                     'vsmin',2.8                 ,... % min sed velocity, km/s
                     'vsstd',0.0                 );  % std of sed velocity for perturbation, km/s

modl.crust = struct(...
    ... thickness of the crust
                     'hmax',60                   ,... %60 max xtal crust thickness, km
                     'hmin',25                   ,... %10 min xtal crust thickness, km
                     'hstd',2.5                    ,... % std of xtal crust thickness, for perturbation, km
                ... gaussian prior probability for crust thickness - mean=30, std=10
                     'h_pprior',@(h) 1,...exp(-(h-30).^2/4.^2),... % prior probability 
    ... shear velocity in the crust
                     'vsmax',4.3                 ,...4.5 % max crust spline velocity, km/s
                     'vsmin',2.5                 ,...3.3 % min crust spline velocity, km/s
                     'vsstd',0.08                ,... % std of crust spline velocity for perturbation, km/s
    ... Vp/Vs in the crust
                     'vpvsmax',1.9               ,...1.9 % min crust vpvs ratio
                     'vpvsmin',1.6               ,...1.65 % max crust vpvs ratio
                     'vpvsstd',0.01              ,... % std of crust vpvs ratio for perturbation, km/s
                ... gaussian prior probability for VpVs - mean=1.8, std=0.05
                     'vpvs_pprior',@(vpvs) 1,...exp(-(vpvs-1.7).^2/0.03.^2),... % prior probability 
    ... Xi in the crust
                     'ximax',1.1                 ,...1.05 % min crust Vs radial anis value
                     'ximin',0.9                 ,...1.00 % min crust Vs radial anis value
                     'xistd',0.01                ,... % std of crust Vs radial anis value
    ... knots in the crust
                     'kdstd',2                   ,... % std of knot movement, for perturbation, km
                     'kmax',7                    ,... % max number of spline knots in crust (inc ends)
                     'kmin',2                    );  % min number of spline knots in crust (inc ends)

modl.mantle = struct(...
    ... shear velocity in the mantle
                     'vsmax',5.1                 ,...4.9 % max mantle spline velocity, km/s
                     'vsmin',3.7                 ,...3.7 % min mantle spline velocity, km/s
                     'vsstd',0.08                ,... % std of mantle spline velocity for perturbation, km/s
    ... Xi in the mantle
                     'ximax',1.0                 ,...1.05 % min mantle Vs radial anis value
                     'ximin',1.0                 ,...1.00 % min mantle Vs radial anis value
                     'xistd',0                   ,... % std of mantle Vs radial anis value
    ... knots in the mantle
                     'kdstd',4                   ,... % std of knot movement, for perturbation, km
                     'kmax',15                   ,... % max number of spline knots in mantle (inc ends)
                     'kmin',5                    );  % max number of spline knots in mantle (inc ends)

modl.data = struct('prior_sigma',struct(                 ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',0.3             ,... %    default
                           'lo',0.2              ,... %    low-f
                           'cms',0.3)            ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',0.2             ,... %    default
                           'lo',0.1))            ,... %    low-f
                  	 'RF',struct(                 ... %  Receiver functions
                    	'Ps',struct(              ... %   P-s data
                           'def',0.2             ,... %    default
                           'lo',0.1              ,... %    low-f
                           'cms',0.3)            ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',0.2             ,... %    default
                           'ccp',0.1             ,... %    ccp stack
                           'lo',0.1))            ,... %    low-f
                  	 'HKstack',struct(            ... %  H-K stack
                    	   'P',0.5)              ,... %    P combination
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)           ,... %    group velocities
                    	'HV',struct(             ... %   Rayleigh wave ellipticity
                           'HVr',0.06)           ,... %    HV ratio
                    	'Lov',struct(             ... %   Love waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)))         ,... %    group velocities
                                                  ...  
                  'min_sigma',struct(             ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',1e-2            ,... %    default
                           'lo',5e-3             ,... %    low-f
                           'cms',1e-3)           ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',1e-3            ,... %    default
                           'lo',1e-3))           ,... %    low-f
                  	 'RF',struct(                 ... %  Receiver functions
                    	'Ps',struct(              ... %   P-s data
                           'def',1e-2            ,... %    default
                           'lo',1e-2             ,... %    low-f
                           'cms',1e-2)           ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',1e-2            ,... %    default
                           'ccp',1e-2             ,... %    ccp stack
                           'lo',1e-2))           ,... %    low-f
                  	 'HKstack',struct(            ... %  H-K stack
                    	   'P',0.2)             ,... %    P combination
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',1e-4            ,... %    phase velocities
                           'grV',1e-4)           ,... %    group velocities
                    	'HV',struct(             ... %   Rayleigh wave ellipticity
                           'HVr',1e-4)           ,... %    HV ratio
                    	'Lov',struct(             ... %   Love waves
                           'phV',1e-4            ,... %    phase velocities
                           'grV',1e-4)))         ,... %    group velocities
                                                  ...  
                  'logstd_sigma',0.05            );   % LOGSTD

                 
%% Forward calc. parms
forc = struct(      'mindV',0.05                 ,... % min delta Vs for layerising
                    'nsamps',2^11                ,... % number of samples (more means longer run time)
                    'PSVorZR','PSV'             ,... % whether to rotate data into PSV or keep in ZR
                    'synthperiod',2              );  % period for propmat response
                
%% Data processing parms
datprocess=struct( 'normdata',true               ,... % normalise data in processing
                   'decdata',false               ,... % decimate data in processing
                   'clipmain',false              ,... % whether to clip the main phase with a taper
                   'clipdaughter',true           ,... % whether to clip the daughter component at the timing of the main phase (to account for improper rotation) == CHEAT
                   'PSV',true                    ,... % PSV (if tru) or ZR
                   'Ps',struct(                   ... % P-s data
                      'Twin'                     ,... %   time window    
                      struct('def',[-2 10],       ... %     default	 
                             'cms',[-2 25])      ,... %     crustal multiples 	  
                      'filtf'                    ,... %   filter frequencies    
                      struct('def',[1e3 1e-3]    ,... %     default	[fhi flo] 
                             'lo',[1/3  1e-3]))  ,... %     low-f [fhi flo]
                   'Sp',struct(                   ... % P-s data
                      'Twin'                     ,... %   time window    
                      struct('def',[-35 4])      ,... %     default 	  
                      'filtf'                    ,... %   filter frequencies    
                      struct('def',[1e3 1e-3]    ,... %     default	[fhi flo] 
                             'lo',[1/5  1e-3]))  ,...
                   'CCP',struct(                  ... % CCP data from Emily
                      'parent_zw',30             ,... %   depth-width of "parent" pulse for CCP    
                      'rayp_S',12.303            ,... %   S wave ray parameter (s/deg) for gcarc=0, edep=0 in iasp91  
                      'surv_Vp_vs',[6.1 3.55]    ,... %   [VP, VS] surface velocity values for P-SV rotation
                      'taperz',10               ,... %   taper width at the edges of the Zwin
                      'Zwin'                     ,... %   depth window    
                      struct('def',[20 250]))     ); %     default	 
%                     
         
%% Model Conditions
cond = struct(  'pos_moho',         true         ,... % No negative moho jumps
                'pos_sed2basement', true         ,... % No negative sed bottom jumps
                'pos_seddV',        true         ,... % Monotonic increase of V(p/s) in sediments
                'pos_crustdV',      false         ,... % Monotonic increase of V(p/s) in crust
                'nobigdVmoh',       true         ,... % No Vs moho jumps exceeding 30%
                'no_moho_sharpgrad',true         ,... % No strong or negative gradients on either side of the Moho
                'noVSgt49',         true         );  % No VS exceeding 4.9 km/s
            
%% Synthetic data parms 
synth = struct( 'gcarcs',[70]                 ,... % average gcarc
                'samprate',10                    ,... % sample rate.
                'noisetype','gaussian'           ,... % noisetype - 'gaussian' or 'real' (in which case drawn from station speficied
                'noise_sigma_SW_Ray',0.015        ,... %0.03 std for random added noise for SWs
                'noise_sigma_SW_Lov',0.015        ,... %0.03 std for random added noise for SWs
                'noise_sigma_SW_HV',0.005        ,... %0.03 std for random added noise for SWs
                'noise_sigma_BW_Sp',0.009        ,... %0.02 std for random added noise for SpRFs
                'noise_sigma_BW_Ps',0.012        ,... %0.02 std for random added noise for PsRFs
                'noise_sigma_RF_Sp',0.009        ,... %0.02 std for random added noise for SpRFs
                'noise_sigma_RF_Ps',0.012        ,... %0.02 std for random added noise for PsRFs
                'surf_Vp_Vs',[6.1 3.55]          ,... % [VP, VS] surface velocity values - if empty, uses True vals
                'SW_Ray_phV_periods',logspace(log10(6),log10(167),22)',...  % Rayleigh wave phV periods
                'SW_Ray_grV_periods',logspace(log10(6),log10(40),10)',...  % Rayleigh wave phV periods
                'SW_Lov_phV_periods',logspace(log10(6),log10(40),10)',...  % Love wave phV periods
                'SW_HV_periods',logspace(log10(8),log10(90),11)',...  % Rayleigh wave HV periods
                'synthperiod',2                  ,...  % period for propmat synth
                'nsamps',[]                      );  % number of samples. 
            
%---RF parameters---%
RFparms = struct([]);
% 'IDRF' for iterative time domain method, or 'ETMTM' for Extended time multitaper method
RFparms(1).method = 'IDRF';
%input: time domain deconvolution
RFparms.gauss_t = 0.5;
RFparms.accept_mis=1e-10; %accepted misfit
RFparms.itmax=50; %number of iterations
%input: freq domain deconvolution
RFparms.TB = 1.5;%1.5; %period*bandwidth
RFparms.NT = 2;%2; %number of tapers
RFparms.tag = 'synth'; %data or synth (synthetic) - Zach uses 'synth' & this gives better results
RFparms.Poverlap = 0.95; %fractional window overlap
%=================================%
synth.RFparms = RFparms;


%% Results parms
res = struct('zatdep',[5:5:modl.maxz]');
            
%% Bundle together
par = struct('inv',inv,'mod',modl,'conditions',cond,'forc',forc,'synth',synth,'datprocess',datprocess,'res',res);
