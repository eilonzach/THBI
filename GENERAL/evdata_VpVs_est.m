function [Vp_surf,Vs_surf] = evdata_VpVs_est(eqar,phases,psv_filtfs,psv_win,evinds_do,SNRminmax,SNR_filtfs,SNR_win,verbose)
%  [Vp_surf,Vs_surf] = evdata_VpVs_est(eqar,phases,filtfs,win,evinds_do,verbose)
% clear all
% close all

fprintf('Calculating station surface Vp and Vs from P-SV transform\n')
if nargin < 1 || isempty(eqar)
    datadir = 'DATA/'; % need final slash
    station = 'RSSD';
    network = 'IU';
    gc_lims_all = [30,75]; % 
    load(sprintf('%sdat_%s_%s_%.0fto%.0f',datadir,station,network,gc_lims_all(1),gc_lims_all(2)));
end

if nargin < 2 || isempty(phases)
    phases = {'Ps','Sp'};
end

if nargin < 3 || isempty(psv_filtfs)
    psv_filtfs = [[0.5;1.5],[0.1;0.75]]; % p,s
end

if nargin < 4 || isempty(psv_win)
    psv_win = [[-5;10],[-10;5]]; % p,s
end

if nargin < 5 || isempty(evinds_do)
    evinds = {[1:eqar.norids]',[1:eqar.norids]'};
else
    evinds = {evinds_do,evinds_do};
end

if nargin < 6 || isempty(SNRminmax)
    psv_SNRmin = 3;
    psv_SNRmax = 20;
else
    psv_SNRmin = SNRminmax(1);
    psv_SNRmax = SNRminmax(2);
end

if nargin < 7 || isempty(SNR_filtfs)
    SNR_filtfs = [0.02;100];
end

if nargin < 8 || isempty(SNR_win)
    SNR_win = [[-85 -15];[-5 5]];   
end

if nargin < 9 || isempty(verbose)
    verbose = True;
end

npoles = 2;
tapertime = 5; % s at beginnning and end of window to taper over
% addpath('matguts');



% prep output structures
[~,~,~,vp_range,vs_range]=surf_vp_vs_est([1,1],0,1,[0 0] );
Vp_est_all = nan(eqar.norids,2);
Vs_est_all = nan(eqar.norids,2);
SNR_est = zeros(eqar.norids,2);
Esurf = zeros(length(vp_range),length(vs_range),eqar.norids,2);


for ip = 1:length(phases)
    fprintf('\n %s \n',phases{ip})
    switch phases{ip}(1)
        case 'P', chp = 1; chd = 2; chidx = 1; tdir = 1; % main phase is Z, align on Z
        case 'S', chp = 2; chd = 1; chidx = 2; tdir =-1; % main phase is R, align on T
    end
    samprate = unique(round(1./diff(eqar.tt(:,:,ip))));
    % ====== flip vertical so it is positive upwards ===== %
    eqar.dataZRT(:,:,1,ip) = -eqar.dataZRT(:,:,1,ip);

    %% SNRest
    fprintf(' > estimating SNR on parent\n')
    % data clean for SNRest
    cp = struct('samprate',samprate,'fhi',SNR_filtfs(2),'flo',SNR_filtfs(1),...
                'pretime',-eqar.tt(1,1,ip),'prex',-SNR_win(1,1)-tapertime,'postx',SNR_win(2,2)+tapertime,...
                'npoles',npoles,'norm',1);
    cp.taperx = tapertime./(cp.prex+cp.postx);
    [ SNRdat,~,~,~,~,SNRtt] = data_clean(  eqar.dataZRT(:,:,chp,ip),cp ); 
    
    % SNR after the filter
    nswin = find(SNRtt>=SNR_win(1,1) & SNRtt < SNR_win(1,2)); % noise window 
    dwin  = find(SNRtt>=SNR_win(2,1) & SNRtt < SNR_win(2,2)); % generous first arrival window 
    SNR_est(:,ip) = max(abs(SNRdat(dwin,:)))./2./rms(SNRdat(nswin,:));
    
    % kill low SNRs
    evinds{ip}(SNR_est(evinds{ip},ip)<psv_SNRmin) = [];
    % cap high SNRs
    SNR_est(SNR_est(:,ip)>psv_SNRmax,ip) = psv_SNRmax;
    
    %gc lims 
%     evinds{ip}(eqar.gcarcs(evinds{ip})>max(gc_lims) | eqar.gcarcs(evinds{ip})<min(gc_lims)) = [] ;

    
    %% window for the Vp and Vs estimate
    win = psv_win(:,ip);
    cp = struct('samprate',samprate,'fhi',psv_filtfs(2,ip),'flo',psv_filtfs(1,ip),...
                'pretime',-eqar.tt(1,1,ip),'prex',-win(1),'postx',win(2),...
                'taperx',1./diff(win),'npoles',npoles,'norm',1);
  
    [ datZ,~,~,~,~,~]  = data_clean(eqar.dataZRT(:,evinds{ip},1,ip),cp ); 
    [ datR,~,~,~,~,tt] = data_clean(eqar.dataZRT(:,evinds{ip},2,ip),cp ); 
    
    % just plot particle motion
%     if ip==1;figure(1);plot(datZ,datR); end

    % For each trace, find the max amplitude in the xcor_window; call
    % this the arrival;
    fprintf(' > estimating surf Vp,Vs\n')
    for ie = 1:length(evinds{ip})
        evind = evinds{ip}(ie);
        [ Vp_est_all(evind,ip),Vs_est_all(evind,ip),Esurf(:,:,evind,ip)] = ...
            surf_vp_vs_est([datZ(:,ie),datR(:,ie)],tt,rayp_sdeg2skm(eqar.rayps(evind,ip)),psv_win(:,ip),0 );
    end
end

%% stack
Estk = zeros(length(vp_range),length(vs_range),3);
for ip=1:2
    for ie = 1:length(evinds{ip})
        evind = evinds{ip}(ie);
        Estk(:,:,ip) = Estk(:,:,ip) + SNR_est(evind,ip)*Esurf(:,:,evind,ip)/mingrid(Esurf(:,:,evind,ip));
    end
    Estk(:,:,ip) = Estk(:,:,ip)./mingrid(Estk(:,:,ip));
end
% relative weight of Ps and Sp:
wts = [median(SNR_est(evinds{1},1)),median(SNR_est(evinds{2},2))];
% stack weighted Ps and Sp E grids
Estk(:,:,3) = wts(1)*Estk(:,:,1) + wts(2)*Estk(:,:,2);
Estk(:,:,3) = Estk(:,:,3)./mingrid(Estk(:,:,3));

[~,x1,y1] = mingrid(Estk(:,:,1));
[~,x2,y2] = mingrid(Estk(:,:,2));
[~,x3,y3] = mingrid(Estk(:,:,3));


%% Different estimates
Vs_surf.comp1_median = median(Vs_est_all(evinds{1},1));
Vs_surf.comp1_mean = mean(Vs_est_all(evinds{1},1));
Vs_surf.comp1_meanwtd = mean_wtd(Vs_est_all(evinds{1},1),SNR_est(evinds{1},1));
Vs_surf.comp1_wtstk = vs_range(x1);
Vp_surf.comp1_median = median(Vp_est_all(evinds{1},1));
Vp_surf.comp1_mean = mean(Vp_est_all(evinds{1},1));
Vp_surf.comp1_meanwtd = mean_wtd(Vp_est_all(evinds{1},1),SNR_est(evinds{1},1));
Vp_surf.comp1_wtstk = vp_range(y1);

Vs_surf.comp2_median = median(Vs_est_all(evinds{2},2));
Vs_surf.comp2_mean = mean(Vs_est_all(evinds{2},2));
Vs_surf.comp2_meanwtd = mean_wtd(Vs_est_all(evinds{2},2),SNR_est(evinds{2},2));
Vs_surf.comp2_wtstk = vs_range(x2);

Vp_surf.comp2_median = median(Vp_est_all(evinds{2},2));
Vp_surf.comp2_mean = mean(Vp_est_all(evinds{2},2));
Vp_surf.comp2_meanwtd = mean_wtd(Vp_est_all(evinds{2},2),SNR_est(evinds{2},2));
Vp_surf.comp2_wtstk = vp_range(y2);

Vs_surf.comp3_wtstk = vs_range(x3);
Vp_surf.comp3_wtstk = vp_range(y3);

% if verbose
% median, mean, wtdmean, stack
fprintf('        median  mean   wtdmean  stack\n')
fprintf('Vs-Ps   %4.2f    %4.2f    %4.2f    %4.2f\n',Vs_surf.comp1_median,Vs_surf.comp1_mean,Vs_surf.comp1_meanwtd,Vs_surf.comp1_wtstk)
fprintf('Vs-Sp   %4.2f    %4.2f    %4.2f    %4.2f\n',Vs_surf.comp2_median,Vs_surf.comp2_mean,Vs_surf.comp2_meanwtd,Vs_surf.comp2_wtstk)
fprintf('Vs-both                         %4.2f\n',Vs_surf.comp3_wtstk)
fprintf('Vp-Ps   %4.2f    %4.2f    %4.2f    %4.2f\n',Vp_surf.comp1_median,Vp_surf.comp1_mean,Vp_surf.comp1_meanwtd,Vp_surf.comp1_wtstk)
fprintf('Vp-Sp   %4.2f    %4.2f    %4.2f    %4.2f\n',Vp_surf.comp2_median,Vp_surf.comp2_mean,Vp_surf.comp2_meanwtd,Vp_surf.comp2_wtstk)
fprintf('Vp-both                         %4.2f\n',Vp_surf.comp3_wtstk)
% end


if verbose
figure(12), clf, set(gcf,'pos',[15 370 1426 428])
subplot(1,3,1), hold on, title('Ps error map')
    contourf(vs_range,vp_range,Estk(:,:,1),40,'edgecolor','none'), colorbar
    scatter(Vs_surf.comp1_wtstk,Vp_surf.comp1_wtstk,50,'r','filled')
subplot(1,3,2), hold on, title('Sp error map')
    contourf(vs_range,vp_range,Estk(:,:,2),40,'edgecolor','none'), colorbar
    scatter(Vs_surf.comp2_wtstk,Vp_surf.comp2_wtstk,50,'g','filled')
subplot(1,3,3), hold on, title('Combined error map')
    contourf(vs_range,vp_range,Estk(:,:,3),40,'edgecolor','none'), colorbar
    scatter(Vs_surf.comp1_wtstk,Vp_surf.comp1_wtstk,50,'r','filled')
    scatter(Vs_surf.comp2_wtstk,Vp_surf.comp2_wtstk,50,'g','filled')
    scatter(Vs_surf.comp3_wtstk,Vp_surf.comp3_wtstk,50,'w','filled')



figure(11)
subplot(1,2,1), plot(SNR_est(:,1),Vs_est_all(:,1),'o')
subplot(1,2,2), plot(SNR_est(:,2),Vs_est_all(:,2),'o')
end


end

