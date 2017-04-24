function [phV,grV] = run_mineos(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% [phV,grV] = run_mineos(model,swperiods,ID,ifdelete,ifplot,ifverbose)
% 
% Function to run the MINEOS for a given model and extract the phase
% velocities at a bunch of input periods. If you keep the output files
% around (ifdelete==false) then they can be used to calculate perturbation
% kernels with the complementary run_kernelcalc.m script


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
qfile = [ID,'.q'];
logfile = [ID,'.log'];

% standard inputs, don't get re-written
modefile = 'safekeeping/modefile.200mhz';
qmod= 'safekeeping/qmod';
%% =======================================================================
wd = pwd;
cd('/Users/zeilon/Documents/MATLAB/BayesianJointInv/matlab_to_mineos');
ifanis = any(model.Sanis) || any(model.Panis);

%% write MINEOS executable and input files format
write_cardfile(cardfile,model.z,model.VP,model.VS,model.rho);
% writeMINEOSmodefile(modefile, ) 
writeMINEOSexecfile( execfile,cardfile,modefile,qmod,eigfile,ofile1,qfile,logfile);
system(['chmod u+x ' execfile]);


%% do MINEOS on it
if ifverbose
    fprintf('    > Running MINEOS normal mode summation code. \n    > Will take some time...')
end
[status,cmdout] = system(['./',execfile]);
if ifverbose
     fprintf(' success!\n')
end
%% read modes output
[phV,grV] = readMINEOS_qfile(qfile,swperiods);
phV = phV(:);
grV = grV(:);


%% delete files
if ifdelete
    delete(execfile,cardfile,eigfile,ofile1,qfile);
    if exist(logfile,'file')==2, delete(logfile); end
end
cd(wd);

%% plot
if ifplot
    figure(88), clf; set(gcf,'pos',[331 385 848 613]);
    ax1 = axes; hold on;
    % dispersion curves
    hd(1)=plot(ax1,swperiods,phV,'o-','linewidth',2);
    hd(2)=plot(ax1,swperiods,grV,'o-','linewidth',2);
    hl = legend(ax1,hd,{'Phase (c)','Group (U)'},'location','southeast');
    set(hl,'fontsize',16,'fontweight','bold')
    set(ax1,'fontsize',16)
    xlabel(ax1,'Period (s)','interpreter','latex','fontsize',22)
    ylabel(ax1,'Velocity (km/s)','interpreter','latex','fontsize',22)
end


 