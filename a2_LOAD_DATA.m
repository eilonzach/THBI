function [trudata,par,stadeets] = a2_LOAD_DATA(par,proj,resdir,avardir)
% function to load data (differs depending on which synth or real)

fprintf('LOADING data\n')
global THBIpath TRUEmodel %#ok<NUSED>

if strcmp(proj.name,'SYNTHETICS')
    
    fprintf(' > Creating custom model and synthetic data\n')
    
    % make synth model
    z0_SYNTH_MODEL_splinemoduse(par,0);  %close(95)
    save([resdir,'/trumodel'],'TRUEmodel');

    % make synth data
    [ trudata ] = z1_SYNTH_DATA(par,0); % in ZRT format
    if strcmp(par.synth.noisetype,'real')
        [ trudata,par ] = z2_NOISIFY_SYNTH( trudata, par,par.synth.noise_sta_deets );
    end

    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);

elseif strcmp(proj.name,'LAB_tests')
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
	try stadeets = irisFetch.Stations('station',par.nwk,par.sta,'*','*'); 
    catch, load([proj.rawdatadir,'/stainfo_master.mat']); 
        stadeets = struct('Latitude',stainfo(strcmp({stainfo.StationCode},par.sta)).Latitude,...
                          'Longitude',stainfo(strcmp({stainfo.StationCode},par.sta)).Longitude);
    end
    
    [trudata,zeroDstr] = load_data(avardir,par.sta,par.nwk,par.gc);
    par.sta = [par.sta,zeroDstr];
    
    % distribute data for different processing (e.g. _lo, _cms)
    trudata = duplicate_data_distribute(trudata,par);
    % set prior sigma as the mean of the data uncertainty if available
    par = prior_sigma_distribute(par,trudata);
    
    % get rid of data that wont be used in inversion - NB NEED EXACT DATA MATCH
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)
        if all(~strcmp(trudtypes{idt},par.inv.datatypes)) % no match
            fprintf('WARNING - removing %s data from trudata\n',trudtypes{idt})
            trudata = rmfield(trudata,trudtypes{idt});
        end
    end
    
    % get rid of remaining data that is not in the right cluster
    trudtypes = fieldnames(trudata);
    for idt = 1:length(trudtypes)    
        if par.inv.BWclust~=0 && any(regexp(trudtypes{idt},'BW'))
            fprintf('WARNING - removing %s data not in cluster %.0f\n',trudtypes{idt},par.inv.BWclust)
            trudata.(trudtypes{idt}) = trudata.(trudtypes{idt})([trudata.(trudtypes{idt}).clust]==par.inv.BWclust);
            % if getting rid of that cluster killed the data type, then
            % delete the data fieldname
            if isempty(trudata.(trudtypes{idt}))
                fprintf('THAT WAS ALL THE %s data\n',trudtypes{idt})
                trudata = rmfield(trudata,trudtypes{idt});
            end
        end
    end
    
    % get rid of datatypes from par if not (any longer) in real data
    kill = false(length(par.inv.datatypes),1);
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; pdt = parse_dtype( dtype ); 
        % since we added in all data types above, if any not left now, it
        % is because they don't exist in the real data - should not use
        if ~isfield(trudata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW') 
            kill(idt) = true;
        end
    end
    par.inv.datatypes(kill) = [];

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

if ~exist('stadeets','var')
        stadeets = struct('Latitude',[],'Longitude',[]);
end
if ~exist('sta','var')
        sta = 'SYN';
end 

end

