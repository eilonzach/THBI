function redoprior = par_equiv_prior(par1,par2)
% redoprior = par_equiv_prior(par1,par2)
%  
% Function to check if two parameter structures are identical as far as the
% parameters that govern the prior - i.e. does a prior have to be
% recalculated?

% default
redoprior = false;

% conditions
if ~isequal(par1.conditions,par2.conditions) 
    redoprior = true; return
end

% model bounds 
if ~isequal(par1.mod.maxz,par2.mod.maxz) 
    redoprior = true; return
end
if ~isequal(par1.mod.maxkz,par2.mod.maxkz) 
    redoprior = true; return
end
if ~isequal(par1.mod.dz,par2.mod.dz) 
    redoprior = true; return
end

% seds
if ~isequal(par1.mod.sed,par2.mod.sed) 
    redoprior = true; return
end
     
% crust
mcfns1 = fieldnames(par1.mod.crust);
mcfns2 = fieldnames(par2.mod.crust);
if ~isequal(mcfns1,mcfns2) 
    redoprior = true;
end
% loop over crust fieldnames
for iff = 1:length(mcfns1)
    mcfn = mcfns1{iff};
    try % will work if not function handle
        if par1.mod.crust.(mcfn) ~= par2.mod.crust.(mcfn)
            redoprior = true;
        end
    catch % must be function handle
        var = strtok(mcfn,'_'); % grab type of variable ('h','vpvs')
        % test with 10 random values
        hh = random('unif',par1.mod.crust.([var,'min']),par1.mod.crust.([var,'max']),10,1);
        if par1.mod.crust.(mcfn)(hh) ~= par2.mod.crust.(mcfn)(hh)
            redoprior = true; return
        end
    end       
end

% seds
if ~isequal(par1.mod.mantle,par2.mod.mantle) 
    redoprior = true; return
end

    
end

