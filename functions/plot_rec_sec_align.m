function plot_rec_sec_align(eqar,phase,comp,filtfs,evinds,dcor,ddpol)
% plot_rec_sec_align(eqar,phase,comp,filtfs,evinds,dcor,ddpol)

if nargin < 5 || isempty(evinds)
    evinds = 1:eqar.norids;
end
if nargin < 6 || isempty(dcor)
    dcor = zeros(length(evinds),1);
end
if nargin < 7 
    ddpol = ones(max(evinds),1);
end

ip = strcmp(eqar.phases,phase);
ic = strcmp(eqar.components,comp);

dt = 1./unique(round(1./diff(eqar.tt(:,1,1))));

% % plot at seaz
% figure(45), clf, set(gcf,'pos',[92 265 1124 781]);
% hold on
% for ie = 1:length(evinds)
%     evind = evinds(ie);
%     dat = eqar.dataZRT(:,evind,ic,ip);
%     dat = dat./max(abs(dat));
%     dat = dat./ddpol(evind);
%     datf = filt_quick(dat,filtfs(1),filtfs(2),dt);
%     datf = datf./max(abs(datf));
%     plot(eqar.tt(:,evind,ip)+dcor(ie),eqar.seazs(evind)+1*datf)
% end

% plot in seaz order
[szs,isz] = sort(eqar.seazs(evinds));
SEP = 1;
SCL = 1.4*SEP;
figure(66), clf, set(gcf,'pos',[92 265 1124 781]);
set(gca,'ytick',SEP*[1:length(evinds)],'yticklabel',num2str(szs,'%.0f'))
hold on
for ie = 1:length(evinds)
    evind = evinds(ie);
    dat = eqar.dataZRT(:,evind,ic,ip);
    dat = dat./max(abs(dat));
    dat = dat./ddpol(evind);
    datf = filt_quick(dat,filtfs(1),filtfs(2),dt);
    datf = datf./max(abs(datf));
    plot(eqar.tt(:,evind,ip)+dcor(ie),SEP*isz(ie)+SCL*datf,'linewidth',1.5)
    text(eqar.tt(end,evind,ip)+dcor(ie)+1,SEP*isz(ie),num2str(evind),'fontsize',15,'fontweight','bold')
end
xlim([eqar.tt(1,evind,ip),eqar.tt(end,evind,ip)]);
ylim([1-2*SEP SEP*length(evinds)+2*SEP])
set(gca,'fontsize',15)
title('Record Section','fontsize',18)
xlabel('Time from arrival (s)','fontsize',17)
ylabel('Sta-evt AZ ( ? )','fontsize',17)















end