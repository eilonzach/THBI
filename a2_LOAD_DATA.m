function [trudata] = a2_LOAD_DATA(par,projname,resdir)
% function to load data (differs depending on which synth or real)

fprintf('LOADING data\n')
global projdir THBIpath TRUEmodel

if strcmp(projname,'SYNTHETICS')
    
    fprintf(' > Creating custom model and synthetic data\n')
    
    % make synth model
    z0_SYNTH_MODEL_fig3_splinemod(par,0);  %close(95)
    save([resdir,'/trumodel'],'TRUEmodel');

    % make synth data
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
    if strcmp(par.synth.noisetype,'real')
        [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,par.synth.noise_sta_deets );
    end

    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
    stadeets = struct('Latitude',[],'Longitude',[]);

elseif strcmp(projname,'LAB_tests')
	z0_SYNTH_MODEL_LAB_TEST(par,par.synth.model.zsed,par.synth.model.zmoh,par.synth.model.zlab,par.synth.model.wlab,par.synth.model.flab,1) ;
	save([resdir,'/trumodel'],'TRUEmodel');

	% make synth data
	[ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
	trudata_noiseless = trudata;

	trudata = trudata_noiseless;
	if strcmp(par.synth.noisetype,'real')
	%     [ trudata,par ] = z2_NOISIFY_SYNTH_makestack( trudata, par,noise_sta_deets );
		[ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,par.synth.noise_sta_deets );
    end
    
    trudata = duplicate_data_distribute(trudata,par);
    stadeets = struct('Latitude',[],'Longitude',[]);

else
	try stadeets = irisFetch.Stations('station',nwk,sta,'*','*'); 
    catch, load([proj.rawdatadir,'/stainfo_master.mat']); 
        stadeets = struct('Latitude',stainfo(strcmp({stainfo.StationCode},sta)).Latitude,...
                          'Longitude',stainfo(strcmp({stainfo.StationCode},sta)).Longitude);
    end

%     [~,~,~,TRUEmodel.Z,TRUEmodel.vs,TRUEmodel.vp,TRUEmodel.rho] = RD_1D_Vprofile; close(gcf);
    [trudata,zeroDstr] = load_data(avardir,sta,nwk,gc);
    sta = [sta,zeroDstr];
    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
    
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
        if ~isfield(trudata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW') 
            trudata.(dtype) = trudata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
        end
        % set prior sigma as geometric mean of data sigma
        if isfield(trudata.(dtype),'sigma') && ~isnan(geomean(trudata.(dtype).sigma))
            par.mod.data.prior_sigma.(pdt{1}).(pdt{2}).(pdt{3}) = geomean(trudata.(dtype).sigma);
        end
    end

    % get rid of data that wont be used in inversion - NB NEED EXACT DATA MATCH
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)
        if all(~strcmp(trudtypes{idt},par.inv.datatypes)) % no match
            fprintf('WARNING - removing %s data from trudata\n',trudtypes{idt})
            trudata = rmfield(trudata,trudtypes{idt});
        end
        if par.inv.BWclust~=0 && any(regexp(trudtypes{idt},'BW'))
            fprintf('WARNING - removing %s data not in cluster %.0f\n',trudtypes{idt},par.inv.BWclust)
            trudata.(trudtypes{idt}) = trudata.(trudtypes{idt})([trudata.(trudtypes{idt}).clust]==par.inv.BWclust);
        end
    end

end

% save data
save([resdir,'/trudata_ORIG'],'trudata');

% window, filter data 
for idt = 1:length(par.inv.datatypes)
    dtype = par.inv.datatypes{idt};
    [ trudata ] = predat_process( trudata,dtype,par);
end
% plot_TRU_WAVEFORMS(trudata);
% plot_TRUvsPRE(trudata,trudata);
save([resdir,'/trudata_USE'],'trudata');
end

