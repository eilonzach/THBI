function [ misfits,par ] = back_compatible( misfits_orig,par_orig )
%[ misfits ] = back_compatible( misfits_orig )
%  
% Script to turn old misfit categories into the new arrangement of
% datatypes

%% MISFITS
misfits = misfits_orig;
%% loop through each chain
nchains = length(misfits);
goodchains = true(nchains,1);

for iii = 1:nchains
    
    if nchains>1, mf = misfits_orig{iii};   else mf = misfits_orig; end
    if isempty(mf), continue; end
    
    mf1 = struct('globmaxL',mf.globmaxL,...
                    'lastL',mf.lastL,...
                    'lastlogL',mf.lastL,...
                    'iter',mf.iter,...
                    'chi2sum',mf.chi2,...
                    'Nstored',mf.Nstored,...
                    'logLike','mf.logLike',...
                    'Like','mf.Like',...
                    'chi2',struct([]),...
                    'rms',struct([]),...
                    'E2',struct([]));
	% chi2
	if isfield(mf,'chi2_SW')
        for jj = 1:length( mf.chi2_SW)
            mf1.chi2(jj,1).SW_Ray_phV = mf.chi2_SW(jj); 
        end
    end
	if isfield(mf,'chi2_ps')
        for jj = 1:length( mf.chi2_ps)
            mf1.chi2(jj,1).BW_Ps = mf.chi2_ps(jj); 
        end
    end    
    if isfield(mf,'chi2_sp')
        for jj = 1:length( mf.chi2_sp)
            mf1.chi2(jj,1).BW_Sp = mf.chi2_sp(jj); 
        end
    end 
	% rmns
	if isfield(mf,'rms_SW')
        for jj = 1:length( mf.rms_SW)
            mf1.rms(jj,1).SW_Ray_phV = mf.rms_SW(jj); 
        end
    end
	if isfield(mf,'rms_ps')
        for jj = 1:length( mf.rms_ps)
            mf1.rms(jj,1).BW_Ps = mf.rms_ps(jj); 
        end
    end    
    if isfield(mf,'rms_sp')
        for jj = 1:length( mf.rms_sp)
            mf1.rms(jj,1).BW_Sp = mf.rms_sp(jj); 
        end
    end 
	% norm
	if isfield(mf,'norm_SW')
        for jj = 1:length( mf.norm_SW)
            mf1.E2(jj,1).SW_Ray_phV = mf.norm_SW(jj); 
        end
    end
	if isfield(mf,'norm_ps')
        for jj = 1:length( mf.norm_ps)
            mf1.E2(jj,1).BW_Ps = mf.norm_ps(jj); 
        end
    end    
    if isfield(mf,'norm_sp')
        for jj = 1:length( mf.norm_sp)
            mf1.E2(jj,1).BW_Sp = mf.norm_sp(jj); 
        end
    end     
    
	misfits{iii} = mf1;

end

%% PAR
par = par_orig;
for id = 1:length(par.inv.datatypes)
    if strcmp(par.inv.datatypes{id},'PsRF')
        par.inv.datatypes{id} = 'BW_Ps';
    end
    if strcmp(par.inv.datatypes{id},'SpRF')
        par.inv.datatypes{id} = 'BW_Sp';
    end
    if strcmp(par.inv.datatypes{id},'SW')
        par.inv.datatypes{id} = 'SW_Ray_phV';
    end    
end


end

