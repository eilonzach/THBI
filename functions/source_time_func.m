function [ stf ] = source_time_func( par )
%SOURCE_TIME_FUNC Summary of this function goes here
%   Detailed explanation goes here

onelayin = struct('zlayt',0,'zlayb',par.mod.maxz,...
                  'Vs',4.1,'Vp',1.81*3.9,'rho',3.1,'nlay',1);

[stf_ps,tt_ps] = run_propmat(onelayin,'mrfunc','Ps',par.synth.samprate, 0.1, par.forc.synthperiod,par.forc.nsamps);
[stf_sp,tt_sp] = run_propmat(onelayin,'mrfunc','Sp',par.synth.samprate, 0.1, par.forc.synthperiod,par.forc.nsamps);

for ic = 1:3
    stf_pow(ic) = diff(tt_ps)'*midpts(stf_ps(:,ic))/max(max(stf_ps));
end
stf_pow

figure(1);
subplot(211); plot(tt_ps,stf_ps), set(gca,'xlim',[0 100])
subplot(212); plot(tt_sp,stf_sp), set(gca,'xlim',[0 100])
end

