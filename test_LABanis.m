% Quick test script to make sure the LAB synthetic model is correctly
% including anisotropy in the forward model...
clear all
close all
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
run([THBIpath,'/a0_STARTUP_BAYES']);
addpath([THBIpath,'/LAB_tests/']);
run('LAB_tests/parms/bayes_inv_parms.m')

zsed = 0;
zmoh = 45;
zlab = 130;
wlab = 40;
flab = 0.05;
par.synth.model = struct('zsed',zsed,'zmoh',zmoh,'zlab',zlab,'wlab',wlab,'flab',flab);
    
    
    
% Lifted straigt from a2_LOAD_DATA
global TRUEmodel
z0_SYNTH_MODEL_LAB_TEST(par,par.synth.model.zsed,par.synth.model.zmoh,par.synth.model.zlab,par.synth.model.wlab,par.synth.model.flab,1) ;

% make synth data
[ trudata1 ] = z1_SYNTH_DATA(par,0); % in ZRT format

% now let's set anis to zero...
TRUEmodel.Sanis = zeros(size(TRUEmodel.Sanis));
[ trudata2 ] = z1_SYNTH_DATA(par,0); % in ZRT format

plot_TRUvsPRE_WAVEFORMS(trudata1,trudata2)

