%% Script to make a bunch of figures for a given station

projname = 'US'; % SYNTHETICS or WYOMING, for now
sta = 'ECSD';
nwk = 'US';
gc = [69,59,40,38,32,66]; % will search for gcarcs +/-3 of this value;
datN = 20;
% baz = 315;

ifsavedat = true;

notes = [...
    '1 chain inversion testing what happens at WVOR. \n',...
    'For this run, we are NOT permitting sediments \n',...
    'For this run, we are using Sp and low-pass Ps \n',...
    'phases along with Rayleigh and AN Love waves/n',...
    'and we are extending the model domain to 300 km depth \n',...
        ];

%% ------------------------- START ------------------------- 
global projdir THBIpath TRUEmodel
THBIpath = '/Users/zeilon/Documents/MATLAB/BayesianJointInv';
projdir = [THBIpath,'/',projname,'/'];
avardir = sprintf('STA_inversions/%s_dat%.0f/',sta,datN);

sta = 'EYMN';
