function [ final_predata ] = final_forward_model( final_model,par,data )
%[ final_predata ] = final_forward_model( final_model,par,data )

ID = 'finalmod';

%% ===================  PUT INTO DATA STRUCTURE  ===================
final_predata = data;

%% ===================  LAYERISE PROFILE  ===================
[zlayt,zlayb,Vslay,Vplay,rholay] = ...
    layerise(final_model.Z,final_model.VSbest,par.forc.mindV,0,final_model.VPbest,final_model.rhobest); 
nlay = length(Vslay);
laymodel = struct('zlayt',zlayt,'zlayb',zlayb,'Vs',Vslay,'Vp',Vplay,'rho',rholay,'nlay',nlay);
if any(isnan(laymodel.rho))
    error('NaN densities')
end


%% =============  PS RFs FROM PROPAGATOR MATRIX  ==============

if any(strcmp(par.inv.datatypes,'PsRF'))
    [ unique_rayps_P,irayps_P ] = rayp_vals( [data.PsRF.rayp] );
    for ir = 1:length(unique_rayps_P)rayp = mean([data.PsRF.rayp]);
        rayp = unique_rayps_P(ir);
        samprate = unique([data.PsRF(irayps_P==ir).samprate]);
        P_inc = rayp2inc(rayp,laymodel.Vp(end),6371-laymodel.zlayb(end));
        [predat_ps,tt_ps] = run_propmat(laymodel,ID,'Ps',samprate, P_inc, par.forc.synthperiod,par.forc.nsamps);
        predat_ps = predat_ps(:,[3,1,2]); % in Z,R,T
        tt_ps_Par = mean(tt_ps(predat_ps(:,1)==min(predat_ps(:,1))));% estimate main P-wave arrival time from first big downswing
        tt_ps = tt_ps - tt_ps_Par;
        tt_ps = round_level(tt_ps,0.001);

        Ps_widewind = [-10 50];
        inwind = (tt_ps >= Ps_widewind(1)) & (tt_ps < Ps_widewind(2)); 
        % crop
        predat_ps = predat_ps(inwind,:);
        tt_ps = tt_ps(inwind);  
        % taper
        predat_ps = flat_hanning_win(tt_ps,predat_ps,Ps_widewind(1),Ps_widewind(2),3); % 3s taper
        % normalise to unit energy
        % normf_ps = predat_ps(:,1)'*predat_ps(:,1) + predat_ps(:,2)'*predat_ps(:,2);
        % predat_ps = predat_ps/sqrt(normf_ps);
        % -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_P==ir);
        for iir = 1:length(inds)
            final_predata.PsRF(inds(iir),1).ZRT=predat_ps;
            final_predata.PsRF(inds(iir),1).tt=tt_ps;
            final_predata.PsRF(inds(iir),1).nsamp = length(final_predata.PsRF(inds(iir)).ZRT);
        end
    end
else
    final_predata.PsRF = [];    
end

%% =============  SP RFs FROM PROPAGATOR MATRIX  ==============
if any(strcmp(par.inv.datatypes,'SpRF'))
    [ unique_rayps_S,irayps_S ] = rayp_vals( [data.SpRF.rayp] );
    for ir = 1:length(unique_rayps_S)
        rayp = unique_rayps_S(ir);
        samprate = unique([data.SpRF(irayps_S==ir).samprate]);
        S_inc = rayp2inc(rayp,laymodel.Vs(end),6371-laymodel.zlayb(end));
        [predat_sp,tt_sp] = run_propmat(laymodel,ID,'Sp',samprate, S_inc, par.forc.synthperiod,par.forc.nsamps);
        predat_sp = predat_sp(:,[3,1,2]); % in Z,R,T
        tt_sp_Sar = mean(tt_sp(predat_sp(:,2)==min(predat_sp(:,2)))); % estimate main S-wave arrival time from first big downswing
        tt_sp = tt_sp - tt_sp_Sar;
        tt_sp = round_level(tt_sp,0.001);

        Sp_widewind = [-50 10];
        inwind = (tt_sp >= Sp_widewind(1)) & (tt_sp < Sp_widewind(2)); 
        % crop
        predat_sp = predat_sp(inwind,:);
        tt_sp = tt_sp(inwind);
        % taper
        predat_sp = flat_hanning_win(tt_sp,predat_sp,Sp_widewind(1),Sp_widewind(2),3); % 3s taper
        % normalise to unit energy
        normf_sp = predat_sp(:,1)'*predat_sp(:,1) + predat_sp(:,2)'*predat_sp(:,2);
        predat_sp = predat_sp/sqrt(normf_sp);
% -----------------  PUT INTO DATA STRUCTURE  -----------------
        inds = find(irayps_S==ir);
        for iir = 1:length(inds)
            final_predata.SpRF(inds(iir),1).ZRT=predat_sp;
            final_predata.SpRF(inds(iir),1).tt=tt_sp;
            final_predata.SpRF(inds(iir),1).nsamp = length(final_predata.SpRF(inds(iir)).ZRT);
        end
    end
else
    final_predata.SpRF = [];    
end

%% ===================  CALCULATE PHASE VELOCITIES  ===================
if any(strcmp(par.inv.datatypes,'SW'))
    %% filenames
    execfile = [ID,'.run_mineos'];
    cardfile = [ID,'.model'];
    eigfile = [ID,'.eig'];
    ofile1 = [ID,'.asc1'];
    qfile = [ID,'.q'];
    logfile = [ID,'.log'];

    % standard inputs, don't get re-written
    modefile = 'safekeeping/modefile.200mhz';
    qmod= 'safekeeping/qmod';
    %% =======================================================================
    wd = pwd;
    cd('/Users/zeilon/Documents/MATLAB/BayesianJointInv/matlab_to_mineos');

    %% write MINEOS executable and input files format
    write_cardfile(cardfile,final_model.Z,final_model.VPbest,final_model.VSbest,final_model.rhobest);
    % writeMINEOSmodefile(modefile, ) 
    writeMINEOSexecfile( execfile,cardfile,modefile,qmod,eigfile,ofile1,qfile,logfile);
    system(['chmod u+x ' execfile]);

    %% do MINEOS on it
    fprintf('    > Running MINEOS normal mode summation code. \n    > Will take some time...')
    [status,cmdout] = system(['./',execfile]);
    fprintf(' success!\n')

    %% read modes output
    [phV,grV] = readMINEOS_qfile(qfile,data.SW.periods);
    predat_cph = phV(:);

    %% delete files
    delete(execfile,cardfile,eigfile,ofile1,qfile);
    if exist(logfile,'file')==2, delete(logfile); end
    cd(wd);


else
   predat_cph = [];
end
% -----------------  PUT INTO DATA STRUCTURE  -----------------
% SW
if any(strcmp(par.inv.datatypes,'SW'))
    final_predata.SW.phV = predat_cph;
end

% 
% % Ps
% final_predata.PsRF = final_predata.PsRF(1);
% final_predata.PsRF.ZRT = predat_ps;
% final_predata.PsRF.tt = tt_ps;
% final_predata.PsRF.nsamp = length(predat_ps);
% % Sp
% final_predata.SpRF = final_predata.SpRF(1);
% final_predata.SpRF.ZRT = predat_sp;
% final_predata.SpRF.tt = tt_sp;
% final_predata.SpRF.nsamp = length(predat_sp);
% % SW
% final_predata.SW.phV = predat_cph;


end

