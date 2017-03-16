%% Description
% in this script all the parameters which define the model and inversion
% are set

%% Inversion parms
inv = struct(    'verbose',true                  ,... % option to spit out more information+plots
                 'niter',10000                    ,... % Number of iterations
                 'burnin',2000                   ,... % don't record results before burnin iterations
                 'cooloff',1000                   ,... % # of iterations over which temperature declines as erf
                 'tempmax',4                     ,... % maximum multiple of all standard deviations
                 'bestNmod2keep',2000             ,... % keep only the best N models in each chain, defined here
                 'kerneltolmax',2.5              ,... % kernel max. tolerance - max norm of perturbation before re-calc kernels
                 'kerneltolmin',0.6              ,... % kernel min. tolerance - norm of perturbation that is totally acceptable
                 'nchains',8                     ,... % number of chains to start in parallel
                 'Kweight',1                     ,... % option to weight SW misfit by fraction of kernel in model space
                 'datatypes',{{'SpRF','PsRF','SW','SpRF_lo','PsRF_lo'}});  % any of {{'SpRF','PsRF','SW','SpRF_lo','PsRF_lo'}}

%% Model parms
mod = struct([]);

mod(1).nstas = 1;
mod.maxz = 300;                                      % maximum depth from ref ellipsoid, km
mod.dz = 3;                                          % depth spacing of model, km

mod.sed = struct(    'hmax',0.0000               ,... %5 max sed layer thickness, km
                     'hmin',0.0000               ,... %0 min sed layer thickness, km
                     'hstd',0                    ,... % std of sed layer thickness for perturbation, km
                     'vsmax',3.43                ,... % max sed velocity, km/s
                     'vsmin',3.42                ,... % min sed velocity, km/s
                     'vsstd',0.0                 );  % std of sed velocity for perturbation, km/s

mod.crust = struct(  'hmax',65                   ,... %60 max xtal crust thickness, km
                     'hmin',30                   ,... %10 min xtal crust thickness, km
                     'hstd',3                    ,... % std of xtal crust thickness, for perturbation, km
                     'vsmax',4.1                 ,...4.5 % max crust spline velocity, km/s
                     'vsmin',3.3                 ,...3.3 % min crust spline velocity, km/s
                     'vsstd',0.1                 ,... % std of crust spline velocity for perturbation, km/s
                     'vpvsmax',1.85              ,...1.9 % min crust vpvs ratio
                     'vpvsmin',1.6               ,...1.65 % max crust vpvs ratio
                     'vpvsstd',0.15              ,... % std of crust vpvs ratio for perturbation, km/s
                     'kdstd',2                   ,... % std of knot movement, for perturbation, km
                     'kmax',5                    ,... % max number of spline knots in crust (inc ends)
                     'kmin',2                    );  % min number of spline knots in crust (inc ends)

mod.mantle = struct( 'vsmax',4.9                 ,...4.9 % max mantle spline velocity, km/s
                     'vsmin',3.7                 ,...3.7 % min mantle spline velocity, km/s
                     'vsstd',0.1                 ,... % std of mantle spline velocity for perturbation, km/s
                     'kdstd',4                   ,... % std of knot movement, for perturbation, km
                     'kmax',20                    ,... % max number of spline knots in mantle (inc ends)
                     'kmin',4                    );  % max number of spline knots in mantle (inc ends)

mod.data = struct(  'prior_sigmaPsRF',0.05       ,... % prior PsRF seismogram sigma
                    'prior_sigmaSpRF',0.05       ,... % prior SpRF seismogram sigma
                    'prior_sigmaSW',0.05         ,... % prior SW phvel sigma
                    'prior_sigmaSpRF_lo',0.02    ,... % prior SpRF seismogram sigma for low f data
                    'prior_sigmaPsRF_lo',0.05    ,... % prior PsRF seismogram sigma for low f data
                    ...
                    'min_sigmaPsRF',0.0002       ,... % minimum PsRF sigma
                    'min_sigmaSpRF',0.0002       ,... % minimum SpRF sigma
                    'min_sigmaSW',0.001          ,... % minimum SW phvel sigma
                    'min_sigmaPsRF_lo',0.0002    ,... % minimum PsRF sigma for low f data
                    'min_sigmaSpRF_lo',0.0002    ,... % minimum SpRF sigma for low f data
                    'sigma_logstd',0.05          );   % std of the log_e of the perturbations to sigma

                 
%% Forward calc. parms
forc = struct(      'mindV',0.05                 ,... % min delta Vs for layerising
                    'nsamps',2^11                ,... % number of samples (more means longer run time)
                    'synthperiod',2              );  % period for propmat synth
                
%% Data processing parms
datprocess=struct( 'normdata',true               ,... % normalise data in processing
                   'decdata',false                ,... % decimate data in processing
                   'Twin'                        ,... % time windows structure...
                     struct('SpRF',[-30 5]       ,... % time window about main S arrival 
                            'PsRF',[-5 25]       ,... % time window about main P arrival 
                            'SpRF_lo',[-30 5]    ,... % time window about main S arrival for low f data
                            'PsRF_lo',[-5 25])   ,... % time window about main P arrival for low f data
                   'filtf'                       ,... % filter structure...
                     struct('SpRF',1./[2 100]    ,... % SpRF filter [fhi flo]
                            'PsRF',1./[1 100]    ,... % PsRF filter [fhi flo]
                            'SpRF_lo',1./[5 100] ,... % SpRF_lo filter [fhi flo]
                            'PsRF_lo',1./[3 100]));   % PsRF_lo filter [fhi flo]
                    
         
%% Model Conditions
cond = struct(  'pos_moho',         true         ,... % No negative moho jumps
                'pos_sed2basement', true         ,... % No negative sed bottom jumps
                'nobigdVmoh',       true         ,... % No Vs moho jumps exceeding 30%
                'pos_crustdV',      false        ,... % Monotonic increase of V(p/s) in crust
                'pos_seddV',        true         ,... % Monotonic increase of V(p/s) in sediments
                'noVSgt49',         true         );  % No VS exceeding 4.9 km/s
            
%% Synthetic data parms 
synth = struct( 'inc',20                         ,... % incidence angle
                'samprate',10                    ,... % sample rate.
                'noise_sigmaSW',0.01             ,... %0.03 std for random added noise for SWs
                'noise_sigmaSpRF',0.005           ,... %0.02 std for random added noise for SpRFs
                'noise_sigmaPsRF',0.008           ,... %0.02 std for random added noise for PsRFs
                'SWperiods',logspace(log10(6),log10(200),25)',...  % surface wave periods
                'nsamps',[]                      );  % number of samples. 

            %% Bundle together
par = struct('inv',inv,'mod',mod,'conditions',cond,'forc',forc,'synth',synth,'datprocess',datprocess);
