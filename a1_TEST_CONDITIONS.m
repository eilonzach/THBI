function [ ifpass ] = a1_TEST_CONDITIONS( model,par,ifverbose )
%[ ifpass ] = a1_TEST_CONDITIONS( model,par )
% 
% Function to check model against conditions laid out in the parameter
% file. These conditions comprise constraints on the model space that limit
% the suite of available models in some intelligent way so as to a priori
% restrict the model space. This function houses all the possible
% constraints, each of which is flagged on/off by an input parameter
% variable
% 
% Available constraints: 
%   "pos_moho" = No negative moho jumps
%   "pos_sed2basement" = No negative sed bottom jumps
%   "nobigdV" = No V(p/s) jumps exceeding 30%
%   "pos_crustdV" = Monotonic increase of V(p/s) in crust
%   "pos_seddV" = Monotonic increase of V(p/s) in sediments
%   "noVSgt49" = No VS exceeding 4.9 km/s
% 
% 
%  INPUTS:
%   model - structure containing the model information, velocity, depth,
%           discontinuities, etc.
%   par   - structure containing all the parameters, one of which will be a
%           "conditions" structure that flags various constraints on and
%           off. 
%   ifverbose - flag to output message with condition that failed

if nargin<3 || isempty(ifverbose)
    ifverbose=0;
end

cond = par.conditions;

ifpass = true; %default is to pass

%% NaN model!
if any(isnan(model.VS))
    ifpass = false;
    if ifverbose,fprintf('Failed: NaN model!\n'); end
    return
end

%% No negative moho jumps
if any(strcmp(fieldnames(cond),'pos_moho')) && cond.pos_moho==true
    if diff(model.VS(model.z==model.zmoh))<0
        ifpass = false;
        if ifverbose,fprintf('Failed: negative moho jump\n'); end
        return
    end
end

%% No negative sed bottom jumps
if any(strcmp(fieldnames(cond),'pos_sed2basement')) && cond.pos_sed2basement==true && model.zsed>0
    if diff(model.VS(model.z==model.zsed))<0
        ifpass = false;
        if ifverbose,fprintf('Failed: negative sed bottom jump\n'); end
        return
    end
end

%% No moho V jumps exceeding 30%
if any(strcmp(fieldnames(cond),'nobigdVmoh')) && cond.nobigdVmoh==true
    if model.fdVSmoh > 30
        ifpass = false;
        if ifverbose,fprintf('Failed: moho fractional dVS>30%\n'); end        
        return
    end
    if model.mantmparm.VS_sp(1) - model.crustmparm.VS_sp(end) > 0.9
        ifpass = false;
        if ifverbose,fprintf('Failed: moho absolute dVS>0.8 km/s\n'); end        
        return
    end
end

%% Monotonic increase of V in crust - use the spline knots for this
if any(strcmp(fieldnames(cond),'pos_crustdV')) && cond.pos_crustdV==true
%     if any(diff(model.VS(model.z<model.zmoh))<0) || ...
%        any(diff(model.VP(model.z<model.zmoh))<0)
    if any(diff(model.crustmparm.VS_sp) < 0)
        ifpass = false;
        if ifverbose,fprintf('Failed: crustal velocity non-increasing\n'); end        
        return
    end
end

%% Monotonic increase of V in sediments
if any(strcmp(fieldnames(cond),'pos_seddV')) && cond.pos_seddV==true
    if diff(model.sedmparm.VS) < 0
        ifpass = false;
        if ifverbose,fprintf('Failed: seds velocity non-increasing\n'); end        
        return
    end
end

%% No VS exceeding 4.9 km/s
if any(strcmp(fieldnames(cond),'noVSgt49')) && cond.noVSgt49==true
    if any(model.VS>4.9)
        ifpass = false;
        if ifverbose,fprintf('Failed: VS exceeds 4.9 km/s\n'); end        
        return
    end
end

%% Velocities in each layer within bounds
is0 = 1:find(model.z==model.zsed,1,'first');  if model.zsed==0, is0 = 0; end
ic0 = is0(end)+1:find(model.z==model.zmoh,1,'first');
im0 = ic0(end)+1:model.Nz;
if model.zsed~=0 && any(model.VS(is0)<par.mod.sed.vsmin)
    ifpass = false;
    if ifverbose,fprintf('Failed: Sed VS < %.2f km/s\n',par.mod.sed.vsmin); end        
    return
end
if model.zsed~=0 && any(model.VS(is0)>par.mod.sed.vsmax)  
    ifpass = false;
    if ifverbose,fprintf('Failed: Sed VS > %.2f km/s\n',par.mod.sed.vsmax); end        
    return
end
% indiv. knots
if any(model.crustmparm.VS_sp<par.mod.crust.vsmin) 
    ifpass = false;
    if ifverbose,fprintf('Failed: A Crust VS spline < %.2f km/s\n',par.mod.crust.vsmin); end        
    return
end
if any(model.crustmparm.VS_sp>par.mod.crust.vsmax) 
    ifpass = false;
    if ifverbose,fprintf('Failed: A Crust VS spline > %.2f km/s\n',par.mod.crust.vsmax); end        
    return
end
if any(model.mantmparm.VS_sp<par.mod.mantle.vsmin)
    ifpass = false;
    if ifverbose,fprintf('Failed: A Mantle VS spline < %.2f km/s\n',par.mod.mantle.vsmin); end        
    return
end
if any(model.mantmparm.VS_sp>par.mod.mantle.vsmax) 
    ifpass = false;
    if ifverbose,fprintf('Failed: A Mantle VS spline > %.2f km/s\n',par.mod.mantle.vsmax); end        
    return
end

if any(model.VS(ic0)<par.mod.crust.vsmin) 
    ifpass = false;
    if ifverbose,fprintf('Failed: Crust VS < %.2f km/s\n',par.mod.crust.vsmin); end        
    return
end
if any(model.VS(ic0)>par.mod.crust.vsmax) 
    ifpass = false;
    if ifverbose,fprintf('Failed: Crust VS > %.2f km/s\n',par.mod.crust.vsmax); end        
    return
end
if any(model.VS(im0)<par.mod.mantle.vsmin) 
    ifpass = false;
    if ifverbose,fprintf('Failed: Mantle VS < %.2f km/s\n',par.mod.mantle.vsmin); end        
    return
end
if any(model.VS(im0)>par.mod.mantle.vsmax) 
    ifpass = false;
    if ifverbose,fprintf('Failed: Mantle VS > %.2f km/s\n',par.mod.mantle.vsmax); end        
    return
end


%% No spline knots below the cutoff, except for the basal one
if any(model.mantmparm.knots(model.mantmparm.knots~=par.mod.maxz)>par.mod.maxkz)
    ifpass = false;
    if ifverbose,fprintf('Failed: Non-basal knot below cutoff\n'); end
    return
end




end

