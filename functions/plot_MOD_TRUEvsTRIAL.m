function  plot_MOD_TRUEvsTRIAL( trumod, testmod )
%PLOT_MOD_TRUEVSTRIAL Summary of this function goes here
%   Detailed explanation goes here

figure(85); clf; set(gcf,'pos',[1294 4 622 529])

subplot(121), hold on;
try plot(trumod.VS,trumod.z,'-k','Linewidth',2);catch, plot(trumod.vs,trumod.Z,'-k','Linewidth',2); end
plot(testmod.VS,testmod.z,'-r','Linewidth',2);
set(gca,'ydir','reverse','fontsize',14);
title('Vs','fontsize',20)
xlabel('Vs (km/s)','fontsize',16)
ylabel('Depth (km)','fontsize',16)
xlim([2.8 4.9])
scatter(2.8*ones(testmod.crustmparm.Nkn,1),testmod.crustmparm.knots,70,'ro','filled')
scatter(2.8*ones(testmod.mantmparm.Nkn,1),testmod.mantmparm.knots,70,'ro','filled')

subplot(122), hold on;
try plot(trumod.VP,trumod.z,'-k','Linewidth',2);catch, plot(trumod.vp,trumod.Z,'-k','Linewidth',2); end
plot(testmod.VP,testmod.z,'-r','Linewidth',2);
set(gca,'ydir','reverse','fontsize',14);
title('Vp','fontsize',20)
xlabel('Vp (km/s)','fontsize',16)
ylabel('Depth (km)','fontsize',16)
xlim([5 9])


end

