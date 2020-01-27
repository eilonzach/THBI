function [ trudata,zeroDstr ] = load_data( par )
%[ trudata,cheatstr ] = load_data( par )

%% no phase velocity data for this nwk yet

%% Surface wave data
seismoddir = '/Volumes/data/models_seismic/';
if ~exist(seismoddir,'dir'), seismoddir = regexprep(seismoddir,'~','/Volumes/zeilon'); end 
addpath(seismoddir);

% -------- Rayleigh waves
[Rperiods,RphV]  = Rph_dispcurve_latlon( avar.slat,avar.slon); % grab composite of AN and EQ
% err = error_curve_EQ_latlon( 1./periods,avar.slat,avar.slon,phVerrordir);

if ~isempty(RphV)
    [Rperiods,iT] = sort(Rperiods);
    RphV = RphV(iT);
    SW_Ray_phV = struct('periods',Rperiods,'phV',RphV,'sigma',[]);
else
    SW_Ray_phV=[];
end

%% no love or HV data for this yet



%% collate
trudata = struct('SW_Ray_phV',SW_Ray_phV);


cd(wd);

