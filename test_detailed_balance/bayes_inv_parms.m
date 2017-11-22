%% Description
%  This script sets all parameters that define the model and inversion
%  for the true inversion with the real data

%% Inversion parms
inv = struct(    'verbose',true                 ,... % option to spit out more information+plots
                 'niter',12000                   ,... % Number of iterations
                 'burnin',3000                   ,... % don't record results before burnin iterations
                 'cooloff',1500                  ,... % # of iterations over which temperature declines as erf
                 'tempmax',4                     ,... % maximum multiple of all standard deviations
                 'bestNmod2keep',2000             ,... % keep only the best N models in each chain, defined here       
                 'kerneltolmax',1.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
                 'kerneltolmed',1.0              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'kerneltolmin',0.5              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'maxnkchain',200                ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'nchains',20                     ,... % number of chains to start in parallel
                 'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
                 'datatypes',{{'BW_Sp','BW_Ps_lo','SW_Ray','SW_Lov'}});  % any of {{'SW','BW_x_y' with x='Sp/Ps' and y=' /lo/fl;}}}

%% Model parms
mod = struct([]);

mod(1).nstas = 1;
mod.maxz = 300;                                      % maximum depth in model from ref ellipsoid, km
mod.maxkz = 250;                                     % maximum depth of deepest non-basal knot, km
mod.dz = 2.5;                                        % depth spacing of model, km

mod.sed = struct(    'hmax',0.0                  ,... %5 max sed layer thickness, km
                     'hmin',0.0                  ,... %0 min sed layer thickness, km
                     'hstd',0.0                  ,... % std of sed layer thickness for perturbation, km
                     'vsmax',3.3                 ,... % max sed velocity, km/s
                     'vsmin',2.5                 ,... % min sed velocity, km/s
                     'vsstd',0.00                 );  % std of sed velocity for perturbation, km/s

mod.crust = struct(  'hmax',60                   ,... %60 max xtal crust thickness, km
                     'hmin',25                   ,... %10 min xtal crust thickness, km
                     'hstd',3                    ,... % std of xtal crust thickness, for perturbation, km
                     'vsmax',4.1                 ,...4.5 % max crust spline velocity, km/s
                     'vsmin',2.7                 ,...3.3 % min crust spline velocity, km/s
                     'vsstd',0.08                ,... % std of crust spline velocity for perturbation, km/s
                     'vpvsmax',1.9               ,...1.9 % min crust vpvs ratio
                     'vpvsmin',1.6               ,...1.65 % max crust vpvs ratio
                     'vpvsstd',0.01              ,... % std of crust vpvs ratio for perturbation, km/s
                     'ximax',1.05                ,...1.05 % min crust Vs radial anis value
                     'ximin',1.0                ,...1.00 % min crust Vs radial anis value
                     'xistd',0.005              ,... % std of crust Vs radial anis value
                     'kdstd',2                   ,... % std of knot movement, for perturbation, km
                     'kmax',6                    ,... % max number of spline knots in crust (inc ends)
                     'kmin',6                    );  % min number of spline knots in crust (inc ends)

mod.mantle = struct( 'vsmax',4.9                 ,...4.9 % max mantle spline velocity, km/s
                     'vsmin',3.7                 ,...3.7 % min mantle spline velocity, km/s
                     'vsstd',0.08                ,... % std of mantle spline velocity for perturbation, km/s
                     'ximax',1.06                ,...1.05 % min mantle Vs radial anis value
                     'ximin',1.0                 ,...1.00 % min mantle Vs radial anis value
                     'xistd',0.005               ,... % std of mantle Vs radial anis value
                     'kdstd',4                   ,... % std of knot movement, for perturbation, km
                     'kmax',6                   ,... % max number of spline knots in mantle (inc ends)
                     'kmin',6                    );  % max number of spline knots in mantle (inc ends)

mod.data = struct('prior_sigma',struct(                 ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',0.2             ,... %    default
                           'lo',0.1              ,... %    low-f
                           'cms',0.3)            ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',0.2             ,... %    default
                           'lo',0.1))            ,... %    low-f
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)           ,... %    group velocities
                    	'Lov',struct(             ... %   Love waves
                           'phV',0.05            ,... %    phase velocities
                           'grV',0.06)))         ,... %    group velocities
                                                  ...  
                  'min_sigma',struct(             ... % PRIOR
                  	 'BW',struct(                 ... %  Body waves
                    	'Ps',struct(              ... %   P-s data
                           'def',2e-4            ,... %    default
                           'lo',2e-4             ,... %    low-f
                           'cms',2e-4)           ,... %    crust multiples
                    	'Sp',struct(              ... %   S-p data
                           'def',2e-4            ,... %    default
                           'lo',2e-4))           ,... %    low-f
                  	 'SW',struct(                 ... %  Surface waves
                    	'Ray',struct(             ... %   Rayleigh waves
                           'phV',1e-4            ,... %    phase velocities
                           'grV',1e-4)           ,... %    group velocities
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
                             'lo',[1/5  1e-3]))  );   %     low-f [fhi flo]
                    
         
%% Model Conditions
cond = struct(  'pos_moho',         true         ,... % No negative moho jumps
                'pos_sed2basement', true         ,... % No negative sed bottom jumps
                'nobigdVmoh',       true         ,... % No Vs moho jumps exceeding 30%
                'pos_crustdV',      true        ,... % Monotonic increase of V(p/s) in crust
                'pos_seddV',        true         ,... % Monotonic increase of V(p/s) in sediments
                'noVSgt49',         true         );  % No VS exceeding 4.9 km/s
            

%% Bundle together
par = struct('inv',inv,'mod',mod,'conditions',cond,'forc',forc,'datprocess',datprocess);
