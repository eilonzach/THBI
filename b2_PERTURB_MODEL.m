 function [ model,ptb,p_bd ] = b2_PERTURB_MODEL( model,par,temp )
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

p_bd = 1; %default, for non b/d moves
model0 = model;

ptbopts = {'ptb_Vmod','ptb_disc','ptb_Mmod','ptb_sigdat'};
% opt_rel_P = [3,2,2,1]; % relative probabilities of altering each one
opt_rel_P = [3,2,2,1];
% opt_rel_P = [0,1,0,0];

opt_plim = cumsum(opt_rel_P)/sum(opt_rel_P);

ptb = [];
while isempty(ptb) % use this while loop to re-do perturbations if the type of perturbation doesn't make sense
model=model0; % just make sure we're back to the original model, if a previous perturbation failed. 

optflag = find(opt_plim>=rand(1),1,'first'); % randomly select which type of perturbation to do
switch ptbopts{optflag} % decide what to modify 

    case 'ptb_Vmod' % modify VELOCITY somewhere in current model
        
        lithlay = {'sed','crust','mantle'};
        lay_rel_P = [1,3,5]; % relative probabilities of altering each one
        % if not perturbing seds...
        if par.mod.sed.vsstd==0, lay_rel_P(1)=0; end
        
        lay_plim = cumsum(lay_rel_P)/sum(lay_rel_P);
        layflag = find(lay_plim>=rand(1),1,'first'); % randomly select which layer to perturb

        
        switch lithlay{layflag}
            
            case 'sed'
                
                ind = randi([1,2]); % random index of V to pertub (but not base of seds)
                std = temp.*par.mod.sed.vsstd; % get std of perturbation
                if std==0, continue; end % don't perturb if no perturbation
                V0 = model.sedmparm.VS(ind);
                vma = par.mod.sed.vsmax;
                vmi = par.mod.sed.vsmin;
                ifgd = false;
                while ifgd==false % only do perturbation within the permitted bounds
                dV = random('norm',0,std,1); % calc. random perturbation
                if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                end
                model.sedmparm.VS(ind) = V0 + dV; % insert perturbed val
                if par.inv.verbose, fprintf('    Changed sed VS(%.0f) from %.2f to %.2f\n',ind,V0,V0+dV);end
                ptb = ['sed_VS_',num2str(ind)];
                
            case 'crust'
                
                V_vpvs_xi_kn_ind = randi(7); % select whether to perturb velocity, vpvs ratio, xi value, or knot location
                
                if V_vpvs_xi_kn_ind<=4 % 57% chance we modify VELOCITY
                    ind = randi([1,model.crustmparm.Nkn]); % random index of V to pertub (but not edges!)
                    std = temp.*par.mod.crust.vsstd; % get std of perturbation
                    V0 = model.crustmparm.VS_sp(ind);
                    vma = par.mod.crust.vsmax;
                    vmi = par.mod.crust.vsmin;
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                    end
                    model.crustmparm.VS_sp(ind) = V0 + dV; % insert perturbed val
                    if par.inv.verbose, fprintf('    Changed crustal VS(%.0f) from %.2f to %.2f\n',ind,V0,V0+dV); end
                    ptb = ['crust_VS_',num2str(ind)];
                elseif V_vpvs_xi_kn_ind==5 % 14% chance we modify vpvs ratio
                    std = temp.*par.mod.crust.vpvsstd; % get std of perturbation
                    V0 = model.crustmparm.vpvs;
                    vma = par.mod.crust.vpvsmax;
                    vmi = par.mod.crust.vpvsmin;
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                    end
                    model.crustmparm.vpvs = V0 + dV; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed crustal vpvs from %.2f to %.2f\n',V0,V0+dV); end
                    ptb = 'crust_vpvs';
                elseif V_vpvs_xi_kn_ind==6 % 14% chance we modify xi value
                    std = temp.*par.mod.crust.xistd; % get std of perturbation
                    V0 = model.crustmparm.xi;
                    vma = par.mod.crust.ximax;
                    vmi = par.mod.crust.ximin;
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                    end
                    model.crustmparm.xi = V0 + dV; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed crustal xi from %.2f to %.2f\n',V0,V0+dV); end
                    ptb = 'crust_xi';
                elseif V_vpvs_xi_kn_ind==7 % 14% chance we modify spline knot location
                    if model.crustmparm.Nkn<=2, continue, end % can't change spline location if only two splines, as can't change edges
                    isp2mod = 1+randi(model.crustmparm.Nkn-2);
                    std = temp.*par.mod.crust.kdstd; % get std of perturbation
                    z0 = model.crustmparm.knots(isp2mod);
                    zma = model.crustmparm.knots(end);
                    zmi = model.crustmparm.knots(1);
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dz = random('norm',0,std,1); % calc. random perturbation
                    if (z0+dz < zma) && (z0+dz > zmi), ifgd = true; end
                    end
                    model.crustmparm.knots(isp2mod) = z0 + dz; % insert perturbed val  
                    model.crustmparm.knots = sort(model.crustmparm.knots);
%                     if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot

                    model.crustmparm.fknots = (model.crustmparm.knots-zmi)/(zma-zmi); % insert perturbed frac  
                    model.crustmparm.splines = make_splines(model.crustmparm.knots,par,model.crustmparm.z_sp);
                    if par.inv.verbose, fprintf('    Moved crustal knot %.0f from %.2f to %.2f\n',isp2mod,z0,z0+dz); end
                    ptb = 'crust_knotmv';
                end % vs, vpvs, xi, or knot
                
            case 'mantle'
                
                V_xi_kn_ind = randi(7); % select whether to perturb velocity or xi value or knot location

                if V_xi_kn_ind<=5 % 71% chance we modify VELOCITY
                    ind = randi([1,model.mantmparm.Nsp]); % random index of V to pertub (but not top of mant)
                    std = temp.*par.mod.mantle.vsstd; % get std of perturbation
                    V0 = model.mantmparm.VS_sp(ind);
                    vma = par.mod.mantle.vsmax;
                    vmi = par.mod.mantle.vsmin;
                    ifgd = false;
                    kkk=0;
                    while ifgd==false % only do perturbation within the permitted bounds
                        kkk=kkk+1;
                        if kkk>100
                            error('stuck in a rut')
                        end
                        dV = random('norm',0,std,1); % calc. random perturbation
                        if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                    end
                    model.mantmparm.VS_sp(ind) = V0 + dV; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed mantle VS(%.0f) from %.2f to %.2f\n',ind,V0,V0+dV); end
                    ptb = ['mantle_VS_',num2str(ind)];
                elseif V_xi_kn_ind==6 % 14% chance we modify xi value
                    std = temp.*par.mod.mantle.xistd; % get std of perturbation
                    V0 = model.mantmparm.xi;
                    vma = par.mod.mantle.ximax;
                    vmi = par.mod.mantle.ximin;
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if (V0+dV <= vma) && (V0+dV >= vmi), ifgd = true; end
                    end
                    model.mantmparm.xi = V0 + dV; % insert perturbed val   
                    if par.inv.verbose, fprintf('    Changed crustal xi from %.2f to %.2f\n',V0,V0+dV); end
                    ptb = 'mantle_xi';
                elseif V_xi_kn_ind==7 % 14% chance we modify spline knot location
                    if model.mantmparm.Nkn<=2, continue, end % can't change spline location if only two knots, as can't change edges
                    isp2mod = 1+randi(model.mantmparm.Nkn-2);
                    std = temp.*par.mod.mantle.kdstd; % get std of perturbation
                    z0 = model.mantmparm.knots(isp2mod);
                    zma = par.mod.maxkz;
                    zmi = model.mantmparm.knots(1);
                    ifgd = false;
                    while ifgd==false % only do perturbation within the permitted bounds
                        dz = random('norm',0,std,1); % calc. random perturbation
                        if (z0+dz < zma) && (z0+dz > zmi), ifgd = true; end
                    end
                    model.mantmparm.knots(isp2mod) = z0 + dz; % insert perturbed val  
                    model.mantmparm.knots = sort(model.mantmparm.knots);
%                     if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot                    

                    model.mantmparm.fknots = (model.mantmparm.knots-zmi)/(zma-zmi); % insert perturbed frac  
                    model.mantmparm.splines = make_splines(model.mantmparm.knots,par,model.mantmparm.z_sp);
                    if par.inv.verbose, fprintf('    Moved mantle knot %.0f from %.2f to %.2f\n',isp2mod,z0,z0+dz); end
                    ptb = 'mantle_knotmv';
                end % vs, xi, or knot
                                
        end % switch on which lith layer to modify velocity of
        
    case 'ptb_disc' % modify DISCONTINUITY in current model
        
        disc = {'sed','moh'};
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
                        vt = model.sedmparm.VS(2);
                        vb = model.crustmparm.VS_sp(1);
                        avV = (vt+vb)/2; % mean V of disc
                        dV0 = vb-vt; % V jump at disc
                        
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.sed.vsstd]);
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                            ddV = random('norm',0,std,1); % calc. random perturbation
                            if (dV0+ddV)>=0, ifgd = true; end
                            if avV - (dV0+ddV)/2 > par.mod.sed.vsmax, ifgd = false; end
                            if avV + (dV0+ddV)/2 < par.mod.crust.vsmin, ifgd = false; end
                        end

                        model.sedmparm.VS(2)      = avV - (dV0+ddV)/2;
                        model.crustmparm.VS_sp(1) = avV + (dV0+ddV)/2;
                        
                        % report
                        if par.inv.verbose, fprintf('    Changed sed/crust dV jump from %.2f to %.2f\n',dV0,dV0+ddV); end
                        ptb = 'sed2crust_dV';
                        
                    case 'avV'
                        vt = model.sedmparm.VS(2);
                        vb = model.crustmparm.VS_sp(1);
                        avV0 = (vt+vb)/2; % mean V of disc
                        dV = vb-vt; % V jump at disc 
                        
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.sed.vsstd]);
                        if std==0, continue; end % don't perturb if no perturbation
                        vmi = model.sedmparm.VS(1);
                        vma = model.crustmparm.VS_sp(2);
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                            davV = random('norm',0,std,1); % calc. random perturbation
                            % if seds have to increase in the crust
%                             if (avV0+davV+dV/2 < vma) && (avV0+davV-dV/2 > vmi), ifgd = true; end 
                            if (avV0+davV-dV/2 > vmi), ifgd = true; end
                            if (avV0+davV) - dV/2 > par.mod.sed.vsmax, ifgd = false; end
                            if (avV0+davV) + dV/2 < par.mod.crust.vsmin, ifgd = false; end
                        end
                        
                        model.sedmparm.VS(2)      = (avV0+davV) - dV/2;
                        model.crustmparm.VS_sp(1) = (avV0+davV) + dV/2;
                        
                        
                        %report
                        if par.inv.verbose, fprintf('    Changed sed/crust avV from %.2f to %.2f\n',avV0,avV0+davV); end
                        ptb = 'sed2crust_avV';
                        
                    case 'h'
                
                        % modify the sed depth
                        std = temp.*par.mod.sed.hstd; % get std of perturbation
                        if std==0, continue; end % don't perturb if no perturbation
                        h0 = model.sedmparm.h;
                        hma = par.mod.sed.hmax; 
                        hmi = par.mod.sed.hmin;
                        if hma==hmi, continue; end % don't perturb sed if no perturbation
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                            dh = random('norm',0,std,1); % calc. random perturbation
                            if (h0+dh <= hma) && (h0+dh >= hmi), ifgd = true; end
                        end 
                        if rem(dh,1)==0, dh = dh+1e-5; end
                        
                        model.crustmparm.knots(1) = h0+dh;
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
                        if dh>0
                            zci = model.z(iczt:iczb);
                            vci = model.VS(iczt:iczb);
                        elseif dh<0
                            zci = [model.zsed+dh;model.z(iczt:iczb)];
                            vci = [model.VS(iczt);model.VS(iczt:iczb);];
                        end

                        model.crustmparm.fknots = (model.crustmparm.knots-cminz)/(cmaxz-cminz);
                        model.crustmparm.Nkn = length(model.crustmparm.knots);
                        model.crustmparm.Nsp = model.crustmparm.Nkn+1;
                        [model.crustmparm.splines,model.crustmparm.VS_sp,model.crustmparm.z_sp]...
                            = make_splines(model.crustmparm.knots,par,zci,vci);     
                        
                        %report
                        if par.inv.verbose, fprintf('    Changed sediment depth from %.2f to %.2f\n',h0,cminz); end
                        ptb = 'sed2crust_h';
                
                end % switch on attribute of disc to modify

            case 'moh'
                
                dattrib = {'dV','avV','h'};
                attflag = randi(3); % equal prob of changing any attribute
                
                switch dattrib{attflag}
                    
                    case 'dV'
                        vt = model.crustmparm.VS_sp(end);
                        vb = model.mantmparm.VS_sp(1);
                        avV = (vt+vb)/2; % mean V of disc
                        dV0 = vb-vt; % V jump at disc
                        
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.mantle.vsstd]);
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                            ddV = random('norm',0,std,1); % calc. random perturbation
                            if (dV0+ddV)>=0, ifgd = true; end
                            if avV - (dV0+ddV)/2 > par.mod.crust.vsmax, ifgd = false; end
                            if avV + (dV0+ddV)/2 < par.mod.mantle.vsmin, ifgd = false; end
                        end

                        model.crustmparm.VS_sp(end) = avV - (dV0+ddV)/2;
                        model.mantmparm.VS_sp(1)    = avV + (dV0+ddV)/2;
                        
                        
                        % report 
                        if par.inv.verbose, fprintf('    Changed Moho dV jump from %.2f to %.2f\n',dV0,dV0+ddV); end
                        ptb = 'Moho_dV';
                        
                    case 'avV'
                        vt = model.crustmparm.VS_sp(end);
                        vb = model.mantmparm.VS_sp(1);
                        avV0 = (vt+vb)/2; % mean V of disc
                        dV = vb-vt; % V jump at disc 
                        
                        std = temp.*mean([par.mod.crust.vsstd,par.mod.mantle.vsstd]);
%                         vmi = model.crustmparm.VS_sp(end-1); % enforce increaseing crustal V
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                            davV = random('norm',0,std,1); % calc. random perturbation
%                             if (avV0+davV-dV/2 > vmi), ifgd = true; end
                            ifgd = true;
                            if (avV0+davV) - dV/2 > par.mod.crust.vsmax, ifgd = false; end
                            if (avV0+davV) + dV/2 < par.mod.mantle.vsmin, ifgd = false; end
                        end
                        
                        model.crustmparm.VS_sp(end) = (avV0+davV) - dV/2;
                        model.mantmparm.VS_sp(1)    = (avV0+davV) + dV/2;
                        
                        
                        % report
                        if par.inv.verbose, fprintf('    Changed Moho avV from %.2f to %.2f\n',avV0,avV0+davV); end
                        ptb = 'Moho_avV';
                        
                    case 'h'
                        % modify the moho depth
                        std = temp.*par.mod.crust.hstd; % get std of perturbation
                        h0 = model.crustmparm.h;
                        hma = par.mod.crust.hmax; 
                        hmi = par.mod.crust.hmin;
                        ifgd = false;
                        while ifgd==false % only do perturbation within the permitted bounds
                        dh = random('norm',0,std,1); % calc. random perturbation
                        if (h0+dh <= hma) && (h0+dh >= hmi), ifgd = true; end
                        end 
                        if rem(dh,1)==0, dh = dh+1e-5; end
                        
                        if dh<0 % wary of moving through a crustal knot
                            model.crustmparm.knots(end) = h0 + dh + model.sedmparm.h;
                            model.crustmparm.knots = sort(model.crustmparm.knots);                        
                            model.mantmparm.knots(1) = model.crustmparm.knots(end);
                        elseif dh>0 % wary of moving through mantle knot
                            model.mantmparm.knots(1) = h0 + dh + model.sedmparm.h;
                            model.mantmparm.knots = sort(model.mantmparm.knots);  
                            model.crustmparm.knots(end) = model.mantmparm.knots(1);
                        end
%                         if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
%                         if any(diff(model.crustmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot

                        
                        % -------- modify splines in crust -------- 
                        cminz = model.crustmparm.knots(1);
                        cmaxz = model.crustmparm.knots(end);
                        model.crustmparm.h = cmaxz - cminz;

                        iczt = find(model.z==model.zsed,1,'last'); % here using the old zsed
                        iczb = find(model.z==model.zmoh,1,'first');% here using the old moho
                        if dh>0
                            zci = [model.z(iczt:iczb);model.zmoh+dh];
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
                        if any(model.mantmparm.VS_sp > par.mod.mantle.vsmax) || any(model.mantmparm.VS_sp < par.mod.mantle.vsmin), continue, end % keep in bounds

                        %report
                        if par.inv.verbose, fprintf('    Changed moho depth from %.2f to %.2f\n',model.sedmparm.h + h0,cmaxz); end
                        ptb = 'Moho_h';
                
                end % switch on attribute of disc to modify
                
        end % switch on which discontinuity to modify
            
    case 'ptb_Mmod'  % add or remove layer from current model
        
        lithlay = {'crust','mantle'};
        lay_rel_P = [1,4]; % relative probabilities of altering each one
        
        lay_plim = cumsum(lay_rel_P)/sum(lay_rel_P);
        layflag = find(lay_plim>=rand(1),1,'first'); % randomly select which layer to perturb
        
        add_or_rm = randi(2);
        
        k = model.mantmparm.Nkn + model.crustmparm.Nkn + 2; % number of layers
     
        switch lithlay{layflag}
            
            case 'crust'
                
                iczt = find(model.z==model.zsed,1,'last');
                iczb = find(model.z==model.zmoh,1,'first');
                zci = model.z(iczt:iczb);
                vci = model.VS(iczt:iczb);
                
                k = model.crustmparm.Nkn;
                
                if add_or_rm == 1 % ADDING a layer
                    if model.crustmparm.Nkn+1>par.mod.crust.kmax, continue, end % can't add more than max knots
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
                    ifgd = false;
                    
                    kk=0;badness=false;
                    while ifgd==false % only do perturbation within the permitted bounds
                    dV = random('norm',0,std,1); % calc. random perturbation
                    if length(isp2mod)>1, dV = [dV/2;-dV/2]; end
                    if all(V0+dV <= vmax) && all(V0+dV >= vmin), ifgd = true; end
                    kk=kk+1;
                    if kk>100, badness=true; break; end
                    end
                    if badness, badness=false; continue; end
                    
                    model.crustmparm.VS_sp(isp2mod) = V0 + dV; % insert perturbed val 
                    % prob of the new velocity - see Bodin 2016 eqn D4
                    if length(isp2mod)>1, absdV = diff(dV); else absdV = abs(dV); end
                    % probability of birth
%                     qV2_vnewk = q_v_given_V( 0,dV,std );
%                     DeltaV = vmax-vmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = (1./qV) * k/(k+1);
                    p_bd = k/(k+1);

                    if par.inv.verbose, fprintf('    Added crustal knot at %.2f km\n',znewk); end
                    ptb = 'crust_kn_add';
                    
                elseif add_or_rm == 2 % REMOVING a layer
                    if model.crustmparm.Nkn-1<par.mod.crust.kmin, continue, end % can't have fewer than min knots
                    oldknots = model.crustmparm.knots;
                    isp2rm = 1+randi(model.crustmparm.Nkn-2); % can't remove an edge one
                    V0 = model.crustmparm.VS_sp(isp2rm); zv0 = model.crustmparm.knots(isp2rm);
                    model.crustmparm.knots(isp2rm) = []; % remove selected knot
                    model.crustmparm.fknots(isp2rm) = []; % remove selected knot fraction
                    model.crustmparm.Nkn = model.crustmparm.Nkn-1;
                    model.crustmparm.Nsp = model.crustmparm.Nsp-1;
                    [model.crustmparm.splines,model.crustmparm.VS_sp] ...
                        = make_splines(model.crustmparm.knots,par,zci,vci);
                    % probability of death
%                     V1 = interp1(zci,model.crustmparm.splines*model.crustmparm.VS_sp,zv0);
%                     qV2_vnewk = q_v_given_V( V0,V1,temp.*par.mod.crust.vsstd );
%                     DeltaV = par.mod.crust.vsmax-par.mod.crust.vsmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = qV * k/(k-1);
                    p_bd = k/(k-1);
                    
                    if par.inv.verbose, fprintf('    Removed crustal knot %.0f\n',isp2rm); end
                    ptb = 'crust_kn_rm';
                end
                
            case 'mantle'
                
                imzt = find(model.z==model.zmoh,1,'last');
                zmi = model.z(imzt:end);
                vmi = model.VS(imzt:end);
                
                k = model.mantmparm.Nkn;

                if add_or_rm == 1 % ADDING a layer
                    if model.mantmparm.Nkn+1>par.mod.mantle.kmax, continue, end % can't add more than max knots
                    model.mantmparm.Nkn = model.mantmparm.Nkn+1;
                    model.mantmparm.Nsp = model.mantmparm.Nsp+1;
                    minz = zmi(2);
                    maxz = zmi(end-1);
                    znewk = random('unif',minz,maxz); % depth of new knot
                    model.mantmparm.knots = sort([model.mantmparm.knots;znewk]); % new vector of knots
%                     if any(diff(model.mantmparm.knots)<par.mod.dz), continue, end % can't be too close to existing knot
                    model.mantmparm.fknots = (model.mantmparm.knots-zmi(1))/(zmi(end)-zmi(1));
                    % new splines from old vels
                    [model.mantmparm.splines,model.mantmparm.VS_sp] ...
                        = make_splines(model.mantmparm.knots,par,zmi,vmi);
                    if any(model.mantmparm.VS_sp > par.mod.mantle.vsmax) || any(model.mantmparm.VS_sp < par.mod.mantle.vsmin), continue, end % keep in bounds% identify biggest new spline @ knot depth
                    isp2mod = spline_weight(model.mantmparm.splines,zmi,znewk);                
                    % make new velocities
                    V0 = model.mantmparm.VS_sp(isp2mod);
                    vmax = par.mod.mantle.vsmax;
                    vmin = par.mod.mantle.vsmin;
                    std = temp.*par.mod.mantle.vsstd; % get std of perturbation
                    ifgd = false;
                    kkk = 0;
                    while ifgd==false % only do perturbation within the permitted bounds
                        kkk=kkk+1;
                        if kkk>100
                            error('stuck in a rut')
                        end
                        dV = random('norm',0,std,1); % calc. random perturbation
                        if length(isp2mod)>1, dV = [-sign(diff(V0))*abs(dV)/2;sign(diff(V0))*abs(dV)/2]; end
                        if all(V0+dV <= vmax) && all(V0+dV >= vmin), ifgd = true; end
                    end
                    model.mantmparm.VS_sp(isp2mod) = V0 + dV; % insert perturbed val    
                    % prob of the new velocity - see Bodin 2016 eqn D4
                    if length(isp2mod)>1, absdV = diff(dV); else absdV = abs(dV); end
                    % probability of birth
%                     qV2_vnewk = q_v_given_V( 0,absdV,std );
%                     DeltaV = vmax-vmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = (1./qV) * k/(k+1);
                    p_bd = k/(k+1);
                    
                    if par.inv.verbose, fprintf('    Added mantle knot at %.2f km\n',znewk); end
                    ptb = 'mantle_kn_add';
                    
                elseif add_or_rm == 2 % REMOVING a layer
                    if model.mantmparm.Nkn-1<par.mod.mantle.kmin, continue, end % can't have fewer than min knots
                    isp2rm = 1+randi(model.mantmparm.Nkn-2); % can't remove an edge one
                    V0 = model.mantmparm.VS_sp(isp2rm); zv0 = model.mantmparm.knots(isp2rm);
                    model.mantmparm.knots(isp2rm) = []; % remove selected knot
                    model.mantmparm.fknots(isp2rm) = []; % remove selected knot fraction
                    model.mantmparm.Nkn = model.mantmparm.Nkn-1;
                    model.mantmparm.Nsp = model.mantmparm.Nsp-1;
                    [model.mantmparm.splines,model.mantmparm.VS_sp] ...
                        = make_splines(model.mantmparm.knots,par,zmi,vmi);
                    if any(model.mantmparm.VS_sp > par.mod.mantle.vsmax) || any(model.mantmparm.VS_sp < par.mod.mantle.vsmin), continue, end % keep in bounds
                    % probability of death
%                     V1 = interp1(zmi,model.mantmparm.splines*model.mantmparm.VS_sp,zv0);
%                     qV2_vnewk = q_v_given_V( V0,V1,temp.*par.mod.mantle.vsstd );
%                     DeltaV = par.mod.mantle.vsmax-par.mod.mantle.vsmin;
%                     qV = qV2_vnewk*DeltaV;
%                     p_bd = qV * k/(k-1);
                    p_bd = k/(k-1);

                    if par.inv.verbose, fprintf('    Removed mantle knot %.0f\n',isp2rm); end
                    ptb = 'mantle_kn_rm';
                end      

        end
        
    case 'ptb_sigdat' % change the sigma value for one of the data types
    
    	datatypes = par.inv.datatypes;
    	dtp2mod = randi(length(datatypes)); % even probability of modifying any
        dtype = datatypes{dtp2mod}; % which data type to modify
        Nd = length(model.datahparm.(['sig_',dtype])); % N of data streams of that dtype
        d2mod = randi(Nd); % which data uncertainty (of the streams) to modify
        ptb = ['sig_',datatypes{dtp2mod},'_',num2str(d2mod)];

        lgstd0 = log10(model.datahparm.(['sig_',dtype])(d2mod));
        ifgd = false;
        while ifgd==false % only do perturbation within the permitted bounds
            dlgstd = random('norm',0,temp.*par.mod.data.logstd_sigma,1); % calc. random perturbation
            lgstd1 = lgstd0+dlgstd;
            model.datahparm.(['sig_',dtype])(d2mod) = 10.^lgstd1; % modify that one data stream's error
            if model.datahparm.(['sig_',dtype])(d2mod) >= get_sigma_min_prior( dtype ,par), ifgd=true; end % check not smaller than minimum 
        end
        if par.inv.verbose, fprintf('    Changed %s data stream %.0f error from %.3f to %.3f\n',dtype,d2mod,10.^lgstd0,10.^lgstd1); end
        
%         switch 
%         
%             case 'SW'
%                 lgstd0 = log10(model.datahparm.sigmaSW);
%                 ifgd = false;
%                 while ifgd==false % only do perturbation within the permitted bounds
%                     dlgstd = random('norm',0,temp.*par.mod.data.sigma_logstd,1); % calc. random perturbation
%                     lgstd1 = lgstd0+dlgstd;
%                     model.datahparm.sigmaSW = 10.^lgstd1;
%                     if model.datahparm.sigmaSW >= par.mod.data.min_sigmaSW, ifgd=true; end
%                 end
%                 if par.inv.verbose, fprintf('    Changed SW data error from %.3f to %.3f\n',10.^lgstd0,10.^lgstd1); end
%             case 'SpRF'
%                 lgstd0 = log10(model.datahparm.sigmaSpRF);
%                 ifgd = false;
%                 while ifgd==false % only do perturbation within the permitted bounds
%                     dlgstd = random('norm',0,temp.*par.mod.data.sigma_logstd,1); % calc. random perturbation
%                     lgstd1 = lgstd0+dlgstd;
%                     model.datahparm.sigmaSpRF = 10.^lgstd1;
%                     if model.datahparm.sigmaSpRF >= par.mod.data.min_sigmaSpRF, ifgd=true; end
%                 end
%                 if par.inv.verbose, fprintf('    Changed SpRF data error from %.3f to %.3f\n',10.^lgstd0,10.^lgstd1); end
%             case 'PsRF'
%                 lgstd0 = log10(model.datahparm.sigmaPsRF);
%                 ifgd = false;
%                 while ifgd==false % only do perturbation within the permitted bounds
%                     dlgstd = random('norm',0,temp.*par.mod.data.sigma_logstd,1); % calc. random perturbation
%                     lgstd1 = lgstd0+dlgstd;
%                     model.datahparm.sigmaPsRF = 10.^lgstd1;
%                     if model.datahparm.sigmaPsRF >= par.mod.data.min_sigmaPsRF, ifgd=true; end
%                 end
%                 if par.inv.verbose, fprintf('    Changed PsRF data error from %.3f to %.3f\n',10.^lgstd0,10.^lgstd1); end
%             case 'SpRF_lo'
%                 lgstd0 = log10(model.datahparm.sigmaSpRF_lo);
%                 ifgd = false;
%                 while ifgd==false % only do perturbation within the permitted bounds
%                     dlgstd = random('norm',0,temp.*par.mod.data.sigma_logstd,1); % calc. random perturbation
%                     lgstd1 = lgstd0+dlgstd;
%                     model.datahparm.sigmaSpRF_lo = 10.^lgstd1;
%                     if model.datahparm.sigmaSpRF_lo >= par.mod.data.min_sigmaSpRF_lo, ifgd=true; end
%                 end
%                 if par.inv.verbose, fprintf('    Changed SpRF_lo data error from %.3f to %.3f\n',10.^lgstd0,10.^lgstd1); end
%             case 'PsRF_lo'
%                 lgstd0 = log10(model.datahparm.sigmaPsRF_lo);
%                 ifgd = false;
%                 while ifgd==false % only do perturbation within the permitted bounds
%                     dlgstd = random('norm',0,temp.*par.mod.data.sigma_logstd,1); % calc. random perturbation
%                     lgstd1 = lgstd0+dlgstd;
%                     model.datahparm.sigmaPsRF_lo = 10.^lgstd1;
%                     if model.datahparm.sigmaPsRF_lo >= par.mod.data.min_sigmaPsRF_lo, ifgd=true; end
%                 end
%                 if par.inv.verbose, fprintf('    Changed PsRF_lo data error from %.3f to %.3f\n',10.^lgstd0,10.^lgstd1); end
%         end
                
end

if any(model.mantmparm.VS_sp<par.mod.mantle.vsmin) || any(model.mantmparm.VS_sp>par.mod.mantle.vsmax) 
    fprintf('Mantle VS exceeds prior bounds\n');
end

end
% USE PERTURBED PARAMS TO MAKE NEW 1D PROFILES ETC.
[ model ] = make_mod_from_parms( model,par );


end

