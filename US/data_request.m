clear all

mkdir request_files
cd request_files
addpath('~/Documents/MATLAB/BayesianJointInv/WYOMING/')

% Parameters"
request_details_all = struct(...
        'phases',	{{'P','S'}},... % need double {{ }}
        'gclims',   [30 75],...
        'maglims',  [5.7 7.2],...
        'samprate', 40,...
        'datwind',  [-100 100]);
request_details_allh = struct(...
        'phases',	{{'P','S'}},... % need double {{ }}
        'gclims',   [30 75],...
        'maglims',  [5.7 7.2],...
        'samprate', 40,...
        'datwind',  [-100 100],...
        'chans','HH?');

% request_details_PS1 = struct(...
%         'phases',	{{'P','S'}},... % need double {{ }}
%         'gclims',   [65 75],...
%         'maglims',  [5.7 7.2],...
%         'samprate', 40,...
%         'datwind',  [-100 100]);
% request_details_P2 = struct(...
%         'phases',	{{'P'}},...     % need double {{ }}
%         'gclims',   [45 55],...
%         'maglims',  [5.7 7.2],...
%         'samprate', 40,...
%         'datwind',  [-100 100]);
return
%% Request data    

% stations_use= {'RSSD','WVOR','MOD'};
% nwks_use = {'IU','US','BK'};
% for is = 1:length(stations_use)
%     evdata2_WAVEFORMS_breqfast(s{is},nwks_use{is},true,false,request_details_all);
% end

% RSSD
evdata2_WAVEFORMS_breqfast('RSSD','IU',true,false,request_details_all);
% WVOR
evdata2_WAVEFORMS_breqfast('WVOR','US',true,false,request_details_all);
% MOD
evdata2_WAVEFORMS_breqfast('MOD','BK',true,false,request_details_all);
% HRV
evdata2_WAVEFORMS_breqfast('HRV','IU',true,false,request_details_all);
% HYB
evdata2_WAVEFORMS_breqfast('HYB','G',true,false,request_details_all);
% ECSD - in the craton
evdata2_WAVEFORMS_breqfast('ECSD','US',true,false,request_details_all);
% EYMN - in the Superior craton
evdata2_WAVEFORMS_breqfast('EYMN','US',true,false,request_details_all);
% ULM - in the Superior craton
evdata2_WAVEFORMS_breqfast('ULM','CN',true,false,request_details_all);

% proposal
evdata2_WAVEFORMS_breqfast('ECHS','YO',true,false,request_details_all);
evdata2_WAVEFORMS_breqfast('B02B','YO',true,false,request_details_allh);
evdata2_WAVEFORMS_breqfast('X04','YO',true,false,request_details_allh);
evdata2_WAVEFORMS_breqfast('D04','YO',true,false,request_details_allh);

cd ..

return
%% Process data

% RSSD
evdata2_WAVEFORMS_breqfast('RSSD','IU',false,true,request_details_all);
% WVOR
evdata2_WAVEFORMS_breqfast('WVOR','US',false,true,request_details_all);
% MOD
evdata2_WAVEFORMS_breqfast('MOD','BK',false,true,request_details_all);
% HRV
evdata2_WAVEFORMS_breqfast('HRV','IU',false,true,request_details_all);
% HYB
evdata2_WAVEFORMS_breqfast('HYB','G',false,true,request_details_all);
% ECSD - in the craton
evdata2_WAVEFORMS_breqfast('ECSD','US',false,true,request_details_all);
% EYMN - in the craton
evdata2_WAVEFORMS_breqfast('EYMN','US',false,true,request_details_all);
% ULM - in the Superior craton
evdata2_WAVEFORMS_breqfast('ULM','CN',false,true,request_details_all);



% proposal
evdata2_WAVEFORMS_breqfast('ECHS','YO',false,true,request_details_all);
evdata2_WAVEFORMS_breqfast('B02B','YO',false,true,request_details_allh);
evdata2_WAVEFORMS_breqfast('X04','YO',false,true,request_details_allh);
evdata2_WAVEFORMS_breqfast('D04','YO',false,true,request_details_allh);
