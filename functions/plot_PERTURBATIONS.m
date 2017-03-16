function [ Nptb_gd,Nptb,ptbname ] = plot_PERTURBATIONS( ptb,ifaccepted )
%PLOT_PERTURBATIONS Summary of this function goes here
%   Detailed explanation goes here

Nptb = struct('lithlay',struct('sed_VS',zeros(2,1),...
                               'crust_VS',zeros(4,1),...
                               'mantle_VS',zeros(10,1),...
                               'crust_vpvs',0),...
              'disc',struct('sed2crust',struct('dV',0,'avV',0,'h',0),...
                            'Moho',struct('dV',0,'avV',0,'h',0)),...
              'hypparm',struct('SpRF_sigma',0,'PsRF_sigma',0,'SW_sigma',0));
Nptb_gd = Nptb;
          
%% MAKE NAMES  

fna = fieldnames(Nptb);
for ia = 1:length(fna)
    fnb = fieldnames(Nptb.(fna{ia}));
    for ib = 1:length(fnb)
        if isstruct(Nptb.(fna{ia}).(fnb{ib}))
            fnc = fieldnames(Nptb.(fna{ia}).(fnb{ib}));
            for ic = 1:length(fnc)
                ptbtry = [fnb{ib},'_',fnc{ic}];
                Nptb.(fna{ia}).(fnb{ib}).(fnc{ic}) = sum(strcmp(ptb,ptbtry));
                Nptb_gd.(fna{ia}).(fnb{ib}).(fnc{ic}) = sum(strcmp(ptb(ifaccepted),ptbtry));
            end
        else
            for ic = 1:length(Nptb.(fna{ia}).(fnb{ib}))
                if length(Nptb.(fna{ia}).(fnb{ib})) > 1
                    nstr = ['_',num2str(ic)];
                else
                    nstr = [];
                end
                ptbtry = [fnb{ib},nstr];
                Nptb.(fna{ia}).(fnb{ib})(ic) = sum(strcmp(ptb,ptbtry));
                Nptb_gd.(fna{ia}).(fnb{ib})(ic) = sum(strcmp(ptb(ifaccepted),ptbtry));
            end
        end
    end
end

% kill excess VS from seds
kill = find(Nptb.lithlay.sed_VS==0);
Nptb.lithlay.sed_VS(kill(kill>1)) = [];
% kill excess VS from crust
kill = find(Nptb.lithlay.crust_VS==0);
Nptb.lithlay.crust_VS(kill(kill>1)) = [];
% kill excess VS from mantle
kill = find(Nptb.lithlay.mantle_VS==0);
Nptb.lithlay.mantle_VS(kill(kill>1)) = [];

%% MAKE PLOT
figure(44), clf, set(gcf,'pos',[680 350 560 748]);
hold on
kk = 1;
ptbname  = {};
fna = fieldnames(Nptb);
for ia = 1:length(fna)
    fnb = fieldnames(Nptb.(fna{ia}));
    for ib = 1:length(fnb)
        if isstruct(Nptb.(fna{ia}).(fnb{ib}))
            fnc = fieldnames(Nptb.(fna{ia}).(fnb{ib}));
            for ic = 1:length(fnc)
                ptbname{kk} = regexprep([fnb{ib},'_',fnc{ic}],'_','-');
                A = Nptb.(fna{ia}).(fnb{ib}).(fnc{ic});
                B = Nptb_gd.(fna{ia}).(fnb{ib}).(fnc{ic});
                plot([0 A],kk*[1,1],'c','linewidth',20)
                plot([0 B],kk*[1,1],'b','linewidth',20)
                kk = kk + 1;
            end
        else
            for ic = 1:length(Nptb.(fna{ia}).(fnb{ib}))
                if length(Nptb.(fna{ia}).(fnb{ib})) > 1
                    nstr = ['_',num2str(ic)];
                else
                    nstr = [];
                end
                ptbname{kk} = regexprep([fnb{ib},nstr],'_','-');
                A = Nptb.(fna{ia}).(fnb{ib})(ic);
                B = Nptb_gd.(fna{ia}).(fnb{ib})(ic);
                plot([0 A],kk*[1,1],'c','linewidth',20)
                plot([0 B],kk*[1,1],'b','linewidth',20)
                kk = kk + 1;
            end
        end
    end
end
set(gca,'ytick',[1:kk-1],'yticklabel',ptbname)



