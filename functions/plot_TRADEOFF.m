function plot_TRADEOFF( posterior )
% plot_TRADEOFF( posterior )
% 
% Plot tradeoffs between posterior model parameters

%% PLOT FINAL MODEL
figure(97); clf; set(gcf,'pos',[120 151 900 850])
ax1 = subplot(221); hold on
ax2 = subplot(222); hold on
ax3 = subplot(223); hold on
ax4 = subplot(224); hold on


%% Mantle velocities 1
zm1 = 1;
zm2 = 6;
x = posterior.VSmantle(:,zm1);
y = posterior.VSmantle(:,zm2);
Gx = [ones(length(y),1),y];
Gy = [ones(length(x),1),x];
bx = Gx\x;
by = Gy\y;
yest = by(1) + by(2)*x;
xest = bx(1) + bx(2)*y;
R2y = 1 - sum((y-yest).^2)/sum((y-mean(y)).^2);
R2x = 1 - sum((x-xest).^2)/sum((x-mean(x)).^2);

plot(ax1,x,y,'or');
% plot(ax1,x,yest,'b','linewidth',1.5)
% plot(ax1,xest,y,'c','linewidth',1.5)
set(ax1,'fontsize',16)
xlabel(ax1,sprintf('VS at %.0f km',posterior.zmantle(zm1)),'Fontsize',18)
ylabel(ax1,sprintf('VS at %.0f km',posterior.zmantle(zm2)),'Fontsize',18)
ll = axis(ax1);
if bx(2)<0
    text(ax1,ll(2) - 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','right')
elseif bx(2)>=0
    text(ax1,ll(1) + 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','left')
end

%% Mantle velocities 2
zm1 = 1;
zm2 = 4;
x = posterior.VSmantle(:,zm1);
y = posterior.VSmantle(:,zm2);
Gx = [ones(length(y),1),y];
Gy = [ones(length(x),1),x];
bx = Gx\x;
by = Gy\y;
yest = by(1) + by(2)*x;
xest = bx(1) + bx(2)*y;
R2y = 1 - sum((y-yest).^2)/sum((y-mean(y)).^2);
R2x = 1 - sum((x-xest).^2)/sum((x-mean(x)).^2);

plot(ax2,x,y,'or');
% plot(ax2,x,yest,'b','linewidth',1.5)
% plot(ax2,xest,y,'c','linewidth',1.5)
set(ax2,'fontsize',16)
xlabel(ax2,sprintf('VS at %.0f km',posterior.zmantle(zm1)),'Fontsize',18)
ylabel(ax2,sprintf('VS at %.0f km',posterior.zmantle(zm2)),'Fontsize',18)
ll = axis(ax2);
if bx(2)<0
    text(ax2,ll(2) - 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','right')
elseif bx(2)>=0
    text(ax2,ll(1) + 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','left')
end

%% Moho depth vs velocity jump
x = posterior.zmoh;
y = posterior.fdVSmoh;
Gx = [ones(length(y),1),y];
Gy = [ones(length(x),1),x];
bx = Gx\x;
by = Gy\y;
yest = by(1) + by(2)*x;
xest = bx(1) + bx(2)*y;
R2y = 1 - sum((y-yest).^2)/sum((y-mean(y)).^2);
R2x = 1 - sum((x-xest).^2)/sum((x-mean(x)).^2);

plot(ax3,x,y,'or');
% plot(ax3,x,yest,'b','linewidth',1.5)
% plot(ax3,xest,y,'c','linewidth',1.5)
set(ax3,'fontsize',16)
xlabel(ax3,'Moho depth','Fontsize',18)
ylabel(ax3,'Moho % dVS','Fontsize',18)
ll = axis(ax3);
if bx(2)<0
    text(ax3,ll(2) - 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','right')
elseif bx(2)>=0
    text(ax3,ll(1) + 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','left')
end

%% Seds depth vs velocity jump
x = posterior.zsed;
y = posterior.fdVSsed;
Gx = [ones(length(y),1),y];
Gy = [ones(length(x),1),x];
bx = Gx\x;
by = Gy\y;
yest = by(1) + by(2)*x;
xest = bx(1) + bx(2)*y;
R2y = 1 - sum((y-yest).^2)/sum((y-mean(y)).^2);
R2x = 1 - sum((x-xest).^2)/sum((x-mean(x)).^2);

plot(ax4,x,y,'or');
% plot(ax4,x,yest,'b','linewidth',1.5)
% plot(ax4,xest,y,'c','linewidth',1.5)
set(ax3,'fontsize',16)
xlabel(ax4,'Seds depth','Fontsize',18)
ylabel(ax4,'Seds % dVS','Fontsize',18)
ll = axis(ax4);
if bx(2)<0
    text(ax4,ll(2) - 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','right')
elseif bx(2)>=0
    text(ax4,ll(1) + 0.05*(ll(2)-ll(1)),ll(4) - 0.07*(ll(4)-ll(3)),...
         sprintf('R^2 = %.2f',mean([R2x,R2y])),...
         'Fontsize',18,'HorizontalAlignment','left')
end

save2pdf(97,'parm_tradeoffs','figs');

end

