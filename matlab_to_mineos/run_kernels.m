function [phV_kernels] = run_kernels(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% [phV_kernels] = run_kernels(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% 
% Function to calculate perturbational phase velocity kernels, having
% previously run MINEOS

if nargin < 3 || isempty(ID)
    ID = 'eg';
end
if nargin < 4 || isempty(ifdelete)
    ifdelete = true;
end
if nargin < 5 || isempty(ifplot)
    ifplot = false;
end
if nargin < 6 || isempty(ifverbose)
    ifverbose = true;
end

%% filenames
if ~ischar(ID), ID = num2str(ID);end
execfile = [ID,'.run_mineos'];
cardfile = [ID,'.model'];
eigfile = [ID,'.eig'];
ofile1 = [ID,'.asc1'];
logfile = [ID,'.log'];
execfile_k = [ID,'.run_kernels'];
stripfile = [ID,'.strip'];
tabfile = [ID,'.table'];
qfile = [ID,'.q'];
kernelfile = [ID,'.cvfrechet'];

% standard inputs, don't get re-written
qmod= 'safekeeping/qmod';

%% =======================================================================
wd = pwd;
cd('/Users/zeilon/Documents/MATLAB/BayesianJointInv/matlab_to_mineos');
ifanis = any(model.Sanis) || any(model.Panis);

%% read modes output
[phV,grV] = readMINEOS_qfile(qfile,swperiods);
phV = phV(:);
grV = grV(:);

%% CALCULATE AND READ IN PERTURBATION KERNELS 
%(frechet derivatves of parm perturbation)
    
%% write kernel calc executable
ikernelfiles = writeKERNELCALCexecfile(swperiods,execfile_k,stripfile,eigfile,qmod,tabfile,qfile,kernelfile,ID,logfile);
system(['chmod u+x ' execfile_k]);

%% do the kernel calculating
if ifverbose
    fprintf('    > Calculting kernels from MINEOS output \n    > Will take some time...')
end
[status,cmdout] = system(['./',execfile_k]);
if ifverbose
    fprintf(' success!\n');
end
%% read 
for ip = 1:length(ikernelfiles)
    phV_kernels{ip,1} = readMINEOS_kernelfile(ikernelfiles{ip});
    phV_kernels{ip,1}.period = swperiods(ip);
end


%% delete files
if ifdelete
    delete(execfile,cardfile,eigfile,ofile1,qfile);
	delete(execfile_k,stripfile,tabfile,kernelfile,[tabfile,'_hdr'],[tabfile,'_hdr.branch']);
	for ip = 1:length(ikernelfiles), delete(ikernelfiles{ip}); end
    if exist(logfile,'file')==2, delete(logfile); end
end
cd(wd);

%% plot
if ifplot
    figure(88), clf; set(gcf,'pos',[331 385 1348 713]);
    ax1 = subplot(3,3,[1,4]); cla, hold on;
    ax2 = subplot(3,3,[2,5,8]); cla, hold on;
    ax3 = subplot(3,3,[3,6,9]); cla, hold on;   
    % kernels
    plot_KERNELS( phV_kernels,ifanis,ax2,ax3 )

    % dispersion curves
    hd(1)=plot(ax1,swperiods,phV,'o-','linewidth',2);
    hd(2)=plot(ax1,swperiods,grV,'o-','linewidth',2);
    hl = legend(ax1,hd,{'Phase (c)','Group (U)'},'location','southeast');
    set(hl,'fontsize',16,'fontweight','bold')
    set(ax1,'fontsize',16)
    xlabel(ax1,'Period (s)','interpreter','latex','fontsize',22)
    ylabel(ax1,'Velocity (km/s)','interpreter','latex','fontsize',22)
end


 