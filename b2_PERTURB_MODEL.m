 function [ model,ptb,P_bd ] = b2_PERTURB_MODEL( model,par,temp )
% [ model_out,ptb,p_bd ] = b2_PERTURB_MODEL( model_in, par,temp )
%   
%  function to perturb the model by choosing one of the fields and then
%  choosing one of the model parameters and then perturbing it about its
%  current value, by an amount dictated by the pre-defined standard
%  deviation multiplied by the temperature (which settles down during the
%  inversion, to produce smaller perturbations as iter time goes on).

% Temperature - a multiple of all the standard deviations that decays with
% increasing iteration so that they settle down to smaller perturbations
% temp = (par.inv.tempmax-1)*erfc(2*(ii-1)./par.inv.cooloff) + 1;

% Four options for model perturbation:
%  ptb_Vmod - change velocity somewhere in current model
%  ptb_disc - change discontinuity depth, dV, or mean velocity in current model
%  ptb_Mmod - add or remove splines from the model (transdimensional)
%  ptb_sigdat - change the noise hyperparameters of the data (hierarchical)
% 
% The outputs areL
%  model   - the new model
%  ptb     - a string with information about which model parm was changed
%  p_bd    - the probability of a birth or death, if relevant

P_bd = 1; %default, for non b/d moves
model0 = model;

ptbopts = {'ptb_Vmod','ptb_disc','ptb_Mmod','ptb_sigdat'};
% opt_rel_P = [3,2,2,1]; % relative probabilities of altering each one
opt_rel_P = [3,2,2,1];

opt_plim = cumsum(opt_rel_P)/sum(opt_rel_P);

ptb = [];
while isempty(ptb) % use this while loop to re-do perturbations if the type of perturbation doesn't make sense
model=model0; % just make sure we're back to the original model, if a previous perturbation failed. 

optflag = find(opt_plim>=rand(1),1,'first'); % randomly select which type of perturbation to do
switch ptbopts{optflag} % decide what to modify 

    case 'ptb_Vmod' % modify VELOCITY somewhere in current model
        
        lithlay = {'sed','crust','mantle'};
        lay_rel_P = [1,3,5]; % relative probabilities of altering each one
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        lay_rel_P = [1,3,5]; % relative probabilities of altering each one
        
        if par.mod.sed.vsstd==0, lay_rel_P(1)=0; end % if not perturbing seds...
        
        lay_plim = cumsum(lay_rel_P)/sum(lay_rel_P);
        layflag = find(lay_plim>=rand(1),1,'first'); % randomly select which layer to perturb

        
        switch lithlay{layflag}
            
            case 'sed'
                
                ind = randi([1,2]); % random index of V to pertub (but not base of seds)
                std = temp.*par.mod.sed.vsstd; % get std of perturbation
                if std==0, continue; end % don't perturb if no perturbation
                ptb = ['sed_VS_',num2str(ind)];

                V0 = model.sedmparm.VS(ind);
                V1 = V0 + random('norm',0,std,1); % calc. random perturbation

                model.sedmparm.VS(ind) = V1; % insert perturbed val
                if par.inv.verbose, fprintf('    Changed sed VS(%.0f) from %.2f to %.2f\n',ind,V0,V1);end
                
                % within bounds?
%                 vma = par.mod.sed.vsmax;
%                 vmi = par.mod.sed.vsmin;
%                 if V1>vma || V1<vmi, p_bd = 0; end

            case 'crust'
                
                val2mod = {'Vs','vpvs','xi','kz'}; % select whether to perturb velocity, vpvs ratio, xi value, or knot location
                val2mod_rel_P = [5,2,2,2]; % relative probabilities of altering each one
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 val2mod_rel_P = [5,2,0,2]; % relative probabilities of altering each one

                val2mod_plim = cumsum(val2mod_rel_P)/sum(val2mod_rel_P);
                valflag = find(val2mod_plim>=rand(1),1,'first'); % randomly select which value to perturb
                
                
                switch val2mod{valflag}

                case 'Vs' % 57% chance we modify VELOCITY
                    ind = randi([1,model.crustmparm.Nsp]); % random index of V to pertub (but not edges!)
                    std = temp.*par.mod.crust.vsstd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = ['crust_VS_',num2str(ind)];

                    V0 = model.crustmparm.VS_sp(ind);
                    V1 = V0 + random('norm',0,std,1); % calc. random perturbation
                    
                    model.crustmparm.VS_sp(ind) = V1; % insert perturbed val
                    if par.inv.verbose, fprintf('    Changed crustal VS(%.0f) from %.2f to %.2f\n',ind,V0,V1); end
                    
                    % within bounds?
%                     vma = par.mod.crust.vsmax;
%                     vmi = par.mod.crust.vsmin;
%                     if V1>vma || V1<vmi, p_bd = 0; return, end
                    
                case 'vpvs' % 14% chance we modify vpvs ratio
                    std = temp.*par.mod.crust.vpvsstd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = 'crust_vpvs';
                    
                    V0 = model.crustmparm.vpvs;
                    V1 = V0 + random('norm',0,std,1); % calc. random perturbation

                    model.crustmparm.vpvs = V1; % insert perturbed val
                    if par.inv.verbose, fprintf('    Changed crustal vpvs from %.2f to %.2f\n',V0,V1); end
                    
                    % within bounds?
                    vma = par.mod.crust.vpvsmax;
                    vmi = par.mod.crust.vpvsmin;
                    if V1>vma || V1<vmi, P_bd = 0; return, end
                    
                case 'xi' % 14% chance we modify xi value
                    std = temp.*par.mod.crust.xistd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = 'crust_xi';
                    
                    V0 = model.crustmparm.xi;
                    V1 = V0 + random('norm',0,std,1); % calc. random perturbation
                    
                    model.crustmparm.xi = V1; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed crustal xi from %.2f to %.2f\n',V0,V1); end
                    
                    % within bounds?
                    vma = par.mod.crust.ximax;
                    vmi = par.mod.crust.ximin;
                    if V1>vma || V1<vmi, P_bd = 0; end                    

                case 'kz' % 14% chance we modify spline knot location
                    if model.crustmparm.Nkn<=2, continue, end % can't change spline location if only two splines, as can't change edges
                    isp2mod = 1+randi(model.crustmparm.Nkn-2);
                    std = temp.*par.mod.crust.kdstd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = ['crust_knotmv_',num2str(isp2mod)];
                    
                    z0 = model.crustmparm.knots(isp2mod);
                    z1 = z0 + random('norm',0,std,1); % calc. random perturbation

                    % within bounds?
                    zma = model.crustmparm.knots(end);
                    zmi = model.crustmparm.knots(1);
                    if z1>zma || z1<zmi, P_bd = 0; return, end                    
                    
                    model.crustmparm.knots(isp2mod) = z1; % insert perturbed val  
                    model.crustmparm.knots = sort(model.crustmparm.knots);
%                     if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot

                    model.crustmparm.fknots = (model.crustmparm.knots-zmi)/(zma-zmi); % insert perturbed frac  
                    model.crustmparm.splines = make_splines(model.crustmparm.knots,par,model.crustmparm.z_sp);
                    if par.inv.verbose, fprintf('    Moved crustal knot %.0f from %.2f to %.2f\n',isp2mod,z0,z1); end
                    
                end % vs, vpvs, xi, or knot
                
            case 'mantle'
                
                val2mod = {'Vs','xi','kz'}; % select whether to perturb velocity, xi value, or knot location
                val2mod_rel_P = [5,2,2]; % relative probabilities of altering each one
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 val2mod_rel_P = [5,1,2]; % relative probabilities of altering each one

                val2mod_plim = cumsum(val2mod_rel_P)/sum(val2mod_rel_P);
                valflag = find(val2mod_plim>=rand(1),1,'first'); % randomly select which value to perturb


                switch val2mod{valflag}
                    
                case 'Vs'
                    ind = randi([1,model.mantmparm.Nsp]); % random index of V to pertub (but not top of mant)
                    std = temp.*par.mod.mantle.vsstd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = ['mantle_VS_',num2str(ind)];

                    V0 = model.mantmparm.VS_sp(ind);
                    V1 = V0 + random('norm',0,std,1); % calc. random perturbation
                    
                    model.mantmparm.VS_sp(ind) = V1; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed mantle VS(%.0f) from %.2f to %.2f\n',ind,V0,V1); end
                    
                    % within bounds?
%                     vma = par.mod.mantle.vsmax;
%                     vmi = par.mod.mantle.vsmin;
%                     if V1>vma || V1<vmi, p_bd = 0; end                    
                    
                case 'xi' % 14% chance we modify xi value
                    std = temp.*par.mod.mantle.xistd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = 'mantle_xi';

                    V0 = model.mantmparm.xi;
                    V1 = V0 + random('norm',0,std,1); % calc. random perturbation

                    model.mantmparm.xi = V1; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed crustal xi from %.2f to %.2f\n',V0,V1); end
                    
                    % within bounds?
                    vma = par.mod.mantle.ximax;
                    vmi = par.mod.mantle.ximin;
                    if V1>vma || V1<vmi, P_bd = 0; end                    

                case 'kz' % 14% chance we modify spline knot location
                    if model.mantmparm.Nkn<=2, continue, end % can't change spline location if only two knots, as can't change edges
                    isp2mod = 1+randi(model.mantmparm.Nkn-2);
                    std = temp.*par.mod.mantle.kdstd; % get std of perturbation
                    if std==0, continue; end % don't perturb if no perturbation
                    ptb = ['mantle_knotmv',num2str(isp2mod)];

                    z0 = model.mantmparm.knots(isp2mod);
                    z1 = z0 + random('norm',0,std,1); % calc. random perturbation

                    % within bounds?
                    zma = par.mod.maxkz;
                    zmi = model.mantmparm.knots(1);
                    if z1>zma || z1<zmi, P_bd = 0; return, end                    

                    model.mantmparm.knots(isp2mod) = z1; % insert perturbed val  
                    model.mantmparm.knots = sort(model.mantmparm.knots);
%                     if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot                    

                    model.mantmparm.fknots = (model.mantmparm.knots-zmi)/(zma-zmi); % insert perturbed frac  
                    model.mantmparm.splines = make_splines(model.mantmparm.knots,par,model.mantmparm.z_sp);
                    if par.inv.verbose, fprintf('    Moved mantle knot %.0f from %.2f to %.2f\n',isp2mod,z0,z1); end
                end % vs, xi, or knot
                                
        end % switch on which lith layer to modify velocity of
        
    case 'ptb_disc' % modify DISCONTINUITY in current model
        
        disc = {'sed','moh'};
        disc_rel_P = [1,2]; % relative probabilities of altering each one
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%      %!!!!!!!!!!!!!!!!!!!!!!!!!!!
        disc_rel_P = [1,2]; % relative probabilities of altering each one
        % if not perturbing seds...
        if par.mod.sed.vsstd==0, disc_rel_P(1)=0; end
        
        disc_plim = cumsum(disc_rel_P)/sum(disc_rel_P);
        discflag = find(disc_plim>=rand(1),1,'first'); % randomly select which type of perturbation to do
        
        switch disc{discflag}
        
            case 'sed'
                
                dattrib = {'dV','avV','h'};
                attflag = randi(3); % equal prob of changing any attribute

                switch dattrib{attflag}
                    
                    case 'dV'
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.sed.vsstd]);
                        if std==0, continue; end % don't perturb if no perturbation
                        ptb = 'sed2crust_dV';

                        vt = model.sedmparm.VS(2);
                        vb = model.crustmparm.VS_sp(1);
                        avV = (vt+vb)/2; % mean V of disc
                        dV0 = vb-vt; % V jump at disc
                        dV1 = dV0 + random('norm',0,std,1); % calc. random perturbation
                        
                        if dV1<0, P_bd = 0; return, end
%                         if avV - dV1/2 > par.mod.sed.vsmax, p_bd = 0; return, end
%                         if avV + dV1/2 < par.mod.crust.vsmin, p_bd = 0; return, end

                        model.sedmparm.VS(end)    = avV - dV1/2;
                        model.crustmparm.VS_sp(1) = avV + dV1/2;
                        
                        % report
                        if par.inv.verbose, fprintf('    Changed sed/crust dV jump from %.2f to %.2f\n',dV0,dV1); end

                    case 'avV'
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.sed.vsstd]);
                        if std==0, continue; end % don't perturb if no perturbation
                        ptb = 'sed2crust_avV';

                        vt = model.sedmparm.VS(2);
                        vb = model.crustmparm.VS_sp(1);
                        avV0 = (vt+vb)/2; % mean V of disc
                        dV = vb-vt; % V jump at disc 
                        
                        vsmi = model.sedmparm.VS(1);
                        vcma = model.crustmparm.VS_sp(2);
                        
                        avV1 = avV0 + random('norm',0,std,1); % calc. random perturbation
                        
%                         if (avV1+dV/2 > vcma), p_bd = 0; return, end % if vels have to increase in the crust
                        if (avV1-dV/2 < vsmi), P_bd = 0; return, end
%                         if avV1 - dV/2 > par.mod.sed.vsmax, p_bd = 0; return, end
%                         if avV1 + dV/2 < par.mod.crust.vsmin, p_bd = 0; return, end

                        model.sedmparm.VS(2)      = avV1 - dV/2;
                        model.crustmparm.VS_sp(1) = avV1 + dV/2;
                        
                        %report
                        if par.inv.verbose, fprintf('    Changed sed/crust avV from %.2f to %.2f\n',avV0,avV1); end
                        
                    case 'h'
                
                        % modify the sed depth
                        std = temp.*par.mod.sed.hstd; % get std of perturbation
                        hma = par.mod.sed.hmax; 
                        hmi = par.mod.sed.hmin;
                        if std==0, continue; end % don't perturb if no perturbation
                        if hma==hmi, continue; end % don't perturb sed if no perturbation
                        ptb = 'sed2crust_h';

                        h0 = model.sedmparm.h;
                        
                        h1 = h0 + random('norm',0,std,1); % calc. random perturbation
                        dh = h1-h0;
                        if rem(h1-h0,1)==0, h1 = h1+1e-5; end
                        
                        if h1 < 0 || h1 < hmi, P_bd = 0; return, end
                        if h1 > hma, P_bd = 0; return, end
                        
                        model.crustmparm.knots(1) = h1;
                        model.crustmparm.knots = sort(model.crustmparm.knots);
                        cminz = model.crustmparm.knots(1);
                        cmaxz = model.crustmparm.knots(end);
                        
%                         if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
%                         if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot

                        
                        model.sedmparm.h = cminz;
                        % Need to modify crustal thickness too! Or the moho moves!!
                        model.crustmparm.h = cmaxz-cminz;

                        % modify splines in crust
                        iczt = find(model.z==model.zsed,1,'last');
                        iczb = find(model.z==model.zmoh,1,'first');
                        if (h1-h0)>0
                            zci = model.z(iczt:iczb);
                            vci = model.VS(iczt:iczb);
                        elseif (h1-h0)<0
                            zci = [model.zsed+(h1-h0);model.z(iczt:iczb)];
                            vci = [model.VS(iczt);model.VS(iczt:iczb);];
                        end

                        model.crustmparm.fknots = (model.crustmparm.knots-cminz)/(cmaxz-cminz);
                        model.crustmparm.Nkn = length(model.crustmparm.knots);
                        model.crustmparm.Nsp = model.crustmparm.Nkn+1;
                        [model.crustmparm.splines,model.crustmparm.VS_sp,model.crustmparm.z_sp]...
                            = make_splines(model.crustmparm.knots,par,zci,vci);     
                        % fix any problems in the re-splining
                        model.crustmparm.VS_sp(1) = model0.crustmparm.VS_sp(1);
                        model.crustmparm.VS_sp(model.crustmparm.VS_sp<par.mod.crust.vsmin)=par.mod.crust.vsmin+1e-6; % 1E-6 for rounding errors
                        model.crustmparm.VS_sp(model.crustmparm.VS_sp>par.mod.crust.vsmax)=par.mod.crust.vsmax-1e-6; % 1E-6 for rounding errors
                        
                        %report
                        if par.inv.verbose, fprintf('    Changed sediment depth from %.2f to %.2f\n',h0,cminz); end
                
                end % switch on attribute of disc to modify

            case 'moh'
                
                dattrib = {'dV','avV','h'};
                attflag = randi(3); % equal prob of changing any attribute
                
                switch dattrib{attflag}
                    
                    case 'dV'
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.mantle.vsstd]);
                        if std==0, continue; end % don't perturb if no perturbation
                        ptb = 'Moho_dV';

                        vt = model.crustmparm.VS_sp(end);
                        vb = model.mantmparm.VS_sp(1);
                        avV = (vt+vb)/2; % mean V of disc
                        dV0 = vb-vt; % V jump at disc
                        dV1 = dV0 + random('norm',0,std,1); % calc. random perturbation
                            
                        if dV1<0, P_bd = 0; return, end
%                         if avV - dV1/2 > par.mod.crust.vsmax, p_bd = 0; return, end
%                         if avV + dV1/2 < par.mod.mantle.vsmin, p_bd = 0; return, end
                        
                        model.crustmparm.VS_sp(end) = avV - dV1/2;
                        model.mantmparm.VS_sp(1)    = avV + dV1/2;
                        
                        % report 
                        if par.inv.verbose, fprintf('    Changed Moho dV jump from %.2f to %.2f\n',dV0,dV1); end
                        
                    case 'avV'
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.mantle.vsstd]);
                        if std==0, continue; end % don't perturb if no perturbation
                        vt = model.crustmparm.VS_sp(end);
                        vb = model.mantmparm.VS_sp(1);
                        avV0 = (vt+vb)/2; % mean V of disc
                        dV = vb-vt; % V jump at disc 
                        avV1 = avV0 + random('norm',0,std,1); % calc. random perturbation
                        
%                         vmi = model.crustmparm.VS_sp(end-1); % enforce increaseing crustal V
                        
%                         if avV1 - dV/2 > par.mod.crust.vsmax, p_bd = 0; return, end
%                         if avV1 + dV/2 < par.mod.mantle.vsmin, p_bd = 0; return, end

                        model.crustmparm.VS_sp(end) = avV1 - dV/2;
                        model.mantmparm.VS_sp(1)    = avV1 + dV/2;
                        
                        
                        % report
                        if par.inv.verbose, fprintf('    Changed Moho avV from %.2f to %.2f\n',avV0,avV1); end
                        ptb = 'Moho_avV';
                        
                    case 'h'
                        % modify the moho depth
                        std = temp.*par.mod.crust.hstd; % get std of perturbation
                        if std==0, continue; end % don't perturb if no perturbation
                        ptb = 'Moho_h';

                        h0 = model.crustmparm.h;
                        hma = par.mod.crust.hmax; 
                        hmi = par.mod.crust.hmin;
                        h1 = h0 + random('norm',0,std,1); % calc. random perturbation
                        dh = h1-h0;
                        if rem(dh,1)==0, h1 = h1+1e-5; end
                        
                        if (h1-h0)<0 % wary of moving through a crustal knot
                            model.crustmparm.knots(end) = h1 + model.sedmparm.h;
                            model.crustmparm.knots = sort(model.crustmparm.knots); 
                            h1 = model.crustmparm.knots(end) - model.sedmparm.h;
                            dh = h1-h0;
                            model.mantmparm.knots(1) = model.crustmparm.knots(end);
                        elseif (h1-h0)>0 % wary of moving through mantle knot
                            model.mantmparm.knots(1) = h1 + model.sedmparm.h;
                            model.mantmparm.knots = sort(model.mantmparm.knots); 
                            h1 = model.mantmparm.knots(1) - model.sedmparm.h;
                            dh = h1-h0;
                            model.crustmparm.knots(end) = model.mantmparm.knots(1);
                        end
%                         if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
%                         if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot

                        if h1 < hmi, P_bd = 0; return, end
                        if h1 > hma, P_bd = 0; return, end

                        % -------- modify splines in crust -------- 
                        cminz = model.crustmparm.knots(1);
                        cmaxz = model.crustmparm.knots(end);
                        model.crustmparm.h = cmaxz - cminz;

                        iczt = find(model.z==model.zsed,1,'last'); % here using the old zsed
                        iczb = find(model.z==model.zmoh,1,'first');% here using the old moho
                        if dh>0
                            zci = [model.z(iczt:iczb);model.zmoh+(h1-h0)];
                            vci = [model.VS(iczt:iczb);model.VS(iczb)];
                        elseif dh<0
                            zci = model.z(iczt:iczb);
                            vci = model.VS(iczt:iczb);
                        end
                        
                        model.crustmparm.fknots = (model.crustmparm.knots-cminz)/(cmaxz-cminz);
                        model.crustmparm.Nkn = length(model.crustmparm.knots);
                        model.crustmparm.Nsp = model.crustmparm.Nkn+1;
                        [model.crustmparm.splines,model.crustmparm.VS_sp,model.crustmparm.z_sp]...
                            = make_splines(model.crustmparm.knots,par,zci,vci);    
                        % fix any problems in the re-splining
                        % this comes at the expense of a small misfit in the Vs profile being fit, but is only practically
                        % relevant in the case that the velocity profiles are very close to the prior limits (in which case
                        % maybe expand them)
                        model.crustmparm.VS_sp(end) = model0.crustmparm.VS_sp(end);
                        model.crustmparm.VS_sp(model.crustmparm.VS_sp<par.mod.crust.vsmin)=par.mod.crust.vsmin+1e-6; % 1E-6 for rounding errors
                        model.crustmparm.VS_sp(model.crustmparm.VS_sp>par.mod.crust.vsmax)=par.mod.crust.vsmax-1e-6; % 1E-6 for rounding errors
                        
                        
                        % -------- modify splines in mantle -------- 
                        mminz = model.mantmparm.knots(1);
                        mmaxz = par.mod.maxz + model.selev;

                        imzt = find(model.z==model.zmoh,1,'last');% here using the old moho
                        if dh>0
                            zmi = model.z(imzt:end);
                            vmi = model.VS(imzt:end);
                        elseif dh<0
                             zmi = [cmaxz;model.z(imzt:end)];
                             vmi = [model.VS(imzt);model.VS(imzt:end)];
                        end
                        
                        model.mantmparm.fknots = (model.mantmparm.knots-mminz)/(mmaxz-mminz);
                        model.mantmparm.Nkn = length(model.mantmparm.knots);
                        model.mantmparm.Nsp = model.mantmparm.Nkn+1;
                        [model.mantmparm.splines,model.mantmparm.VS_sp,model.mantmparm.z_sp]...
                            = make_splines(model.mantmparm.knots,par,zmi,vmi);    
                        % fix any problems in the re-splining
                        % this comes at the expense of a small misfit in the Vs profile being fit, but is only practically
                        % relevant in the case that the velocity profiles are very close to the prior limits (in which case
                        % maybe expand them)
                        model.mantmparm.VS_sp(1) = model0.mantmparm.VS_sp(1);
                        model.mantmparm.VS_sp(model.mantmparm.VS_sp<par.mod.mantle.vsmin)=par.mod.mantle.vsmin+1e-6; % 1E-6 for rounding errors
                        model.mantmparm.VS_sp(model.mantmparm.VS_sp>par.mod.mantle.vsmax)=par.mod.mantle.vsmax-1e-6; % 1E-6 for rounding errors

%                         if any(model.mantmparm.VS_sp > par.mod.mantle.vsmax) || any(model.mantmparm.VS_sp < par.mod.mantle.vsmin), p_bd = 0; return, end % keep in bounds

                        %report
                        if par.inv.verbose, fprintf('    Changed moho depth from %.2f to %.2f\n',model.sedmparm.h + h0,cmaxz); end
                
                end % switch on attribute of disc to modify
                
        end % switch on which discontinuity to modify
            
    case 'ptb_Mmod'  % add or remove layer from current model
        
        lithlay = {'crust','mantle'};
        lay_rel_P = [1,4]; % relative probabilities of altering each one
        
        lay_plim = cumsum(lay_rel_P)/sum(lay_rel_P);
        layflag = find(lay_plim>=rand(1),1,'first'); % randomly select which layer to perturb
        
        add_or_rm = randi(2);
        
     
        switch lithlay{layflag}
            
            case 'crust'
                
                iczt = find(model.z==model.zsed,1,'last');
                iczb = find(model.z==model.zmoh,1,'first');
                zci = model.z(iczt:iczb);
                vci = model.VS(iczt:iczb);
                
                k = model.crustmparm.Nkn;
                
                if add_or_rm == 1 % ADDING a layer
                    if par.mod.crust.kmin==par.mod.crust.kmax, continue; end % don't perturb if no perturbation allowed
                    ptb = 'crust_kn_add';
                    if model.crustmparm.Nkn+1>par.mod.crust.kmax, P_bd = 0; return, end % can't add more than max knots
                    % add knot
                    model.crustmparm.Nkn = model.crustmparm.Nkn+1;
                    model.crustmparm.Nsp = model.crustmparm.Nsp+1;
                    % bounds on where new knot is
                    minz = zci(2);
                    maxz = zci(end-1);
                    znewk = random('unif',minz,maxz); % depth of new knot
                    model.crustmparm.knots = sort([model.crustmparm.knots;znewk]); % new vector of knots
%                     if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
                    model.crustmparm.fknots = (model.crustmparm.knots-zci(1))/(zci(end)-zci(1));
                    % new splines from old vels
                    [model.crustmparm.splines,model.crustmparm.VS_sp] ...
                        = make_splines(model.crustmparm.knots,par,zci,vci);
                    % identify biggest new spline @ knot depth
                    isp2mod = spline_weight(model.crustmparm.splines,zci,znewk);
                    % make new velocities
                    V0 = model.crustmparm.VS_sp(isp2mod);
                    vmax = par.mod.crust.vsmax;
                    vmin = par.mod.crust.vsmin;
                    std = temp.*par.mod.crust.vsstd; % get std of perturbation
                    
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if length(isp2mod)>1, dV = [dV/2;-dV/2]; end
                    V1 = V0 + dV;
                    
%                     if any(V1 > vmax), p_bd = 0; return; end
%                     if any(V1 < vmin), p_bd = 0; return; end
                    
                    model.crustmparm.VS_sp(isp2mod) = V1; % insert perturbed val 

                    % prob of the new velocity - see Bodin 2016 eqn D4
                    if length(isp2mod)>1, absdV = diff(dV); else absdV = abs(dV); end
                    % probability of birth
%                     qV2_vnewk = q_v_given_V( 0,dV,std );
%                     DeltaV = vmax-vmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = (1./qV) * k/(k+1);
                    P_bd = k/(k+1);

                    if par.inv.verbose, fprintf('    Added crustal knot at %.2f km\n',znewk); end
                    
                elseif add_or_rm == 2 % REMOVING a layer
                    if par.mod.crust.kmin==par.mod.crust.kmax, continue; end % don't perturb if no perturbation allowed
                    ptb = 'crust_kn_rm';
                    if model.crustmparm.Nkn-1<par.mod.crust.kmin, P_bd = 0; return, end % can't have fewer than min knots
                    oldknots = model.crustmparm.knots;
                    isp2rm = 1+randi(model.crustmparm.Nkn-2); % can't remove an edge one
                    V0 = model.crustmparm.VS_sp(isp2rm); zv0 = model.crustmparm.knots(isp2rm);
                    model.crustmparm.knots(isp2rm) = []; % remove selected knot
                    model.crustmparm.fknots(isp2rm) = []; % remove selected knot fraction
                    model.crustmparm.Nkn = model.crustmparm.Nkn-1;
                    model.crustmparm.Nsp = model.crustmparm.Nsp-1;
                    [model.crustmparm.splines,model.crustmparm.VS_sp] ...
                        = make_splines(model.crustmparm.knots,par,zci,vci);
                    % fix any problems in the re-splining
                    model.crustmparm.VS_sp([1,end]) = model0.crustmparm.VS_sp([1,end]);
                    model.crustmparm.VS_sp(model.crustmparm.VS_sp<par.mod.crust.vsmin)=par.mod.crust.vsmin+1e-6; % 1E-6 for rounding errors
                    model.crustmparm.VS_sp(model.crustmparm.VS_sp>par.mod.crust.vsmax)=par.mod.crust.vsmax-1e-6; % 1E-6 for rounding errors
                    
                    % probability of death
%                     V1 = interp1(zci,model.crustmparm.splines*model.crustmparm.VS_sp,zv0);
%                     qV2_vnewk = q_v_given_V( V0,V1,temp.*par.mod.crust.vsstd );
%                     DeltaV = par.mod.crust.vsmax-par.mod.crust.vsmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = qV * k/(k-1);
                    P_bd = k/(k-1);
                    
                    if par.inv.verbose, fprintf('    Removed crustal knot %.0f\n',isp2rm); end
                end
                
            case 'mantle'
                
                imzt = find(model.z==model.zmoh,1,'last');
                zmi = model.z(imzt:end);
                vmi = model.VS(imzt:end);
                
                k = model.mantmparm.Nkn;

                if add_or_rm == 1 % ADDING a layer
                    if par.mod.mantle.kmin==par.mod.mantle.kmax, continue; end % don't perturb if no perturbation allowed
                    ptb = 'mantle_kn_add';
                    if model.mantmparm.Nkn+1>par.mod.mantle.kmax, P_bd = 0; return, end % can't add more than max knots
                    model.mantmparm.Nkn = model.mantmparm.Nkn+1;
                    model.mantmparm.Nsp = model.mantmparm.Nsp+1;
                    minz = zmi(2);
                    maxz = par.mod.maxkz;
                    znewk = random('unif',minz,maxz); % depth of new knot
                    model.mantmparm.knots = sort([model.mantmparm.knots;znewk]); % new vector of knots
%                     if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
                    model.mantmparm.fknots = (model.mantmparm.knots-zmi(1))/(zmi(end)-zmi(1));
                    % new splines from old vels
                    [model.mantmparm.splines,model.mantmparm.VS_sp] ...
                        = make_splines(model.mantmparm.knots,par,zmi,vmi);

                    isp2mod = spline_weight(model.mantmparm.splines,zmi,znewk);                
                    % make new velocities
                    V0 = model.mantmparm.VS_sp(isp2mod);
                    vmax = par.mod.mantle.vsmax;
                    vmin = par.mod.mantle.vsmin;
                    std = temp.*par.mod.mantle.vsstd; % get std of perturbation
                    
                    dV = random('norm',0,std,1); % calc. random perturbation
%                     if length(isp2mod)>1, dV = [-sign(diff(V0))*abs(dV)/2;sign(diff(V0))*abs(dV)/2]; end
                    if length(isp2mod)>1, dV = [dV/2;-dV/2]; end
                    V1 = V0 + dV;

%                     if any(V1 > vmax), p_bd = 0; return, end
%                     if any(V0 < vmin), p_bd = 0; return, end
                                        
                    model.mantmparm.VS_sp(isp2mod) = V1; % insert perturbed val 

                    % prob of the new velocity - see Bodin 2016 eqn D4
                    if length(isp2mod)>1, absdV = diff(dV); else absdV = abs(dV); end
                    % probability of birth
%                     qV2_vnewk = q_v_given_V( 0,absdV,std );
%                     DeltaV = vmax-vmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = (1./qV) * k/(k+1);
                    P_bd = k/(k+1);
                    
                    if par.inv.verbose, fprintf('    Added mantle knot at %.2f km\n',znewk); end
                    
                elseif add_or_rm == 2 % REMOVING a layer
                    if par.mod.mantle.kmin==par.mod.mantle.kmax, continue; end % don't perturb if no perturbation allowed
                    ptb = 'mantle_kn_rm';
                    if model.mantmparm.Nkn-1<par.mod.mantle.kmin, P_bd = 0; return, end % can't have fewer than min knots
                    isp2rm = 1+randi(model.mantmparm.Nkn-2); % can't remove an edge one
                    V0 = model.mantmparm.VS_sp(isp2rm); zv0 = model.mantmparm.knots(isp2rm);
                    model.mantmparm.knots(isp2rm) = []; % remove selected knot
                    model.mantmparm.fknots(isp2rm) = []; % remove selected knot fraction
                    model.mantmparm.Nkn = model.mantmparm.Nkn-1;
                    model.mantmparm.Nsp = model.mantmparm.Nsp-1;
                    [model.mantmparm.splines,model.mantmparm.VS_sp] ...
                        = make_splines(model.mantmparm.knots,par,zmi,vmi);
                    % fix any problems in the re-splining
                    model.mantmparm.VS_sp([1,end]) = model0.mantmparm.VS_sp([1,end]);
                    model.mantmparm.VS_sp(model.mantmparm.VS_sp<par.mod.mantle.vsmin)=par.mod.mantle.vsmin+1e-6; % 1E-6 for rounding errors
                    model.mantmparm.VS_sp(model.mantmparm.VS_sp>par.mod.mantle.vsmax)=par.mod.mantle.vsmax-1e-6; % 1E-6 for rounding errors

                    % probability of death
%                     V1 = interp1(zmi,model.mantmparm.splines*model.mantmparm.VS_sp,zv0);
%                     qV2_vnewk = q_v_given_V( V0,V1,temp.*par.mod.mantle.vsstd );
%                     DeltaV = par.mod.mantle.vsmax-par.mod.mantle.vsmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = qV * k/(k-1);
                    P_bd = k/(k-1);

                    if par.inv.verbose, fprintf('    Removed mantle knot %.0f\n',isp2rm); end
                end      

        end
        
    case 'ptb_sigdat' % change the sigma value for one of the data types
    
        std = par.mod.data.logstd_sigma;
        if std==0, continue; end % don't perturb if no perturbation allowed
    	datatypes = par.inv.datatypes;
    	dtp2mod = randi(length(datatypes)); % even probability of modifying any
        dtype = datatypes{dtp2mod}; % which data type to modify
        Nd = length(model.datahparm.(['sig_',dtype])); % N of data streams of that dtype
        d2mod = randi(Nd); % which data uncertainty (of the streams) to modify
        ptb = ['sig_',datatypes{dtp2mod},'_',num2str(d2mod)];

        lgstd0 = log10(model.datahparm.(['sig_',dtype])(d2mod));
        
        lgstd1 = lgstd0 + random('norm',0,temp.*std,1); % calc. random perturbation
        
        model.datahparm.(['sig_',dtype])(d2mod) = 10.^lgstd1; % modify that one data stream's error
            
        if model.datahparm.(['sig_',dtype])(d2mod) < get_sigma_min_prior( dtype,par), P_bd = 0; return, end % check not smaller than minimum 

        if par.inv.verbose, fprintf('    Changed %s data stream %.0f error from %.3f to %.3f\n',dtype,d2mod,10.^lgstd0,10.^lgstd1); end
              
end

% if any(model.mantmparm.VS_sp<par.mod.mantle.vsmin) || any(model.mantmparm.VS_sp>par.mod.mantle.vsmax) 
%     p_bd = 0;
% end

end

% no need to recalc model if only hyperparm change or if model will reject
if strcmp(ptbopts{optflag},'ptb_sigdat'), return; end
if P_bd == 0, return; end


%% USE PERTURBED PARAMS TO MAKE NEW 1D PROFILES ETC.
[ model ] = make_mod_from_parms( model,par );


% if strcmp(ptb,'Moho_h') && a1_TEST_CONDITIONS( model, par ) == 0
%     figure(2); clf,hold on
%     plot(model0.z(1:40),model0.VS(1:40),'k',model.z(1:40),model.VS(1:40),'r')
%     plot(model2.z(1:40),model2.VS(1:40),'g')
%     plot(zci,vci,'b',zmi,vmi,'c','linewidth',2)
%     xlim([0 70])
% end

end

