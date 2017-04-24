function [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt )
% [ misfit ] = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt)
%  
% Calculate cross-convolution missfit between observed data and
% (unit-normalised) predicted data

if nargin < 4 || isempty(ifplot)
    ifplot = 0;
end
if nargin < 5 || isempty(SWwt)
    SWwt = ones(size(trudata.SW.phV));
end

%% Ps misfit
if any(strcmp(par.inv.datatypes,'PsRF')) 
    misfit2_ps = zeros(length(trudata.PsRF),1);stfpow_ps_tru=misfit2_ps;stfpow_ps_pre=misfit2_ps;
    for itr = 1:length(trudata.PsRF)
        [ misfit2_ps(itr) ] = xconv_misfit(trudata.PsRF(itr).ZRT(:,1),...
                                           trudata.PsRF(itr).ZRT(:,2),...
                                           predata.PsRF(itr).ZRT(:,1),...
                                           predata.PsRF(itr).ZRT(:,2));
        stfpow_ps_tru(itr) = norm(trudata.PsRF(itr).ZRT(:,1)) +...
                             norm(trudata.PsRF(itr).ZRT(:,2)) +...
                             norm(trudata.PsRF(itr).ZRT(:,3));
        stfpow_ps_pre(itr) = norm(predata.PsRF(itr).ZRT(:,1)) +...
                             norm(predata.PsRF(itr).ZRT(:,2)) +...
                             norm(predata.PsRF(itr).ZRT(:,3));
    end
    misfit2_ps = misfit2_ps./stfpow_ps_tru./stfpow_ps_pre;   
else
    misfit2_ps = 0;
end
%% Sp misfit
if any(strcmp(par.inv.datatypes,'SpRF'))
    misfit2_sp = zeros(length(trudata.SpRF),1);stfpow_sp_tru=misfit2_sp;stfpow_sp_pre=misfit2_sp;
    for itr = 1:length(trudata.SpRF)
        [ misfit2_sp(itr) ] = xconv_misfit(trudata.SpRF(itr).ZRT(:,1),...
                                           trudata.SpRF(itr).ZRT(:,2),...
                                           predata.SpRF(itr).ZRT(:,1),...
                                           predata.SpRF(itr).ZRT(:,2));
        stfpow_sp_tru(itr) = norm(trudata.SpRF(itr).ZRT(:,1)) +...
                             norm(trudata.SpRF(itr).ZRT(:,2)) +...
                             norm(trudata.SpRF(itr).ZRT(:,3));
        stfpow_sp_pre(itr) = norm(predata.SpRF(itr).ZRT(:,1)) +...
                             norm(predata.SpRF(itr).ZRT(:,2)) +...
                             norm(predata.SpRF(itr).ZRT(:,3));
    end
    misfit2_sp = misfit2_sp./stfpow_sp_tru./stfpow_sp_pre;
else
    misfit2_sp = 0;
end

%% Ps_lo misfit
if any(strcmp(par.inv.datatypes,'PsRF_lo'))
    misfit2_ps_lo = zeros(length(trudata.PsRF_lo),1);stfpow_ps_lo_tru=misfit2_ps_lo;stfpow_ps_lo_pre=misfit2_ps_lo;
    for itr = 1:length(trudata.PsRF_lo)
        [ misfit2_ps_lo(itr) ] = xconv_misfit(trudata.PsRF_lo(itr).ZRT(:,1),...
                                              trudata.PsRF_lo(itr).ZRT(:,2),...
                                              predata.PsRF_lo(itr).ZRT(:,1),...
                                              predata.PsRF_lo(itr).ZRT(:,2));
        stfpow_ps_lo_tru(itr) = norm(trudata.PsRF_lo(itr).ZRT(:,1)) +...
                                norm(trudata.PsRF_lo(itr).ZRT(:,2)) +...
                                norm(trudata.PsRF_lo(itr).ZRT(:,3));
        stfpow_ps_lo_pre(itr) = norm(predata.PsRF_lo(itr).ZRT(:,1)) +...
                                norm(predata.PsRF_lo(itr).ZRT(:,2)) +...
                                norm(predata.PsRF_lo(itr).ZRT(:,3));
    end
    misfit2_ps_lo = misfit2_ps_lo./stfpow_ps_lo_tru./stfpow_ps_lo_pre;
else
    misfit2_ps_lo = 0;
end

%% Sp_lo misfit
if any(strcmp(par.inv.datatypes,'SpRF_lo'))
    misfit2_sp_lo = zeros(length(trudata.SpRF_lo),1);stfpow_sp_lo_tru=misfit2_sp_lo;stfpow_sp_lo_pre=misfit2_sp_lo;
    for itr = 1:length(trudata.SpRF_lo)
        [ misfit2_sp_lo(itr) ] = xconv_misfit(trudata.SpRF_lo(itr).ZRT(:,1),...
                                              trudata.SpRF_lo(itr).ZRT(:,2),...
                                              predata.SpRF_lo(itr).ZRT(:,1),...
                                              predata.SpRF_lo(itr).ZRT(:,2));
        stfpow_sp_lo_tru(itr) = norm(trudata.SpRF_lo(itr).ZRT(:,1)) +...
                                norm(trudata.SpRF_lo(itr).ZRT(:,2)) +...
                                norm(trudata.SpRF_lo(itr).ZRT(:,3));
        stfpow_sp_lo_pre(itr) = norm(predata.SpRF_lo(itr).ZRT(:,1)) +...
                                norm(predata.SpRF_lo(itr).ZRT(:,2)) +...
                                norm(predata.SpRF_lo(itr).ZRT(:,3));
    end
    misfit2_sp_lo = misfit2_sp_lo./stfpow_sp_lo_tru./stfpow_sp_lo_pre;
else
    misfit2_sp_lo = 0;
end

%% SW misfit
if any(strcmp(par.inv.datatypes,'SW'))    
    e = (trudata.SW.phV - predata.SW.phV);
    misfit2_SW = e'*(SWwt.*e);
else
    misfit2_SW = 0;
end


%% OVERALL
misfit = struct('PsRF',misfit2_ps,'SpRF',misfit2_sp,'SW',misfit2_SW,'SpRF_lo',misfit2_sp_lo,'PsRF_lo',misfit2_ps_lo);

if ifplot
    plot_TRUvsPRE( trudata.PsRF.ZRT,trudata.PsRF.ZRT,predata.PsRF.ZRT,predata.PsRF.ZRT,par)
end

end
