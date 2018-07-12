function [ parsed_dtype ] = parse_dtype( dtype )
%[ parsed_dtype ] = parse_dtype( dtype )
%   function to read a string that describes a data type and parse it into
%   fields with an order that correspond to the data type. If any fields
%   are un-specified, the value for that field is 'def' (default);
% 
%  if BW, the fields are: { Phase | Twin | filtf }
%           'Phase' can be 'Ps' or 'Sp'
%           'Twin' can be default or 'cms', which indicates crustal multiples
%               included for the P wave
%           'filtf' can be default (no-filter) or 'lo'
%  if SW, the fields are: { RL | PGE }
%           'RL' can be 'Ray' or 'Lov' for Rayleigh or Love waves
%           'PGE' can be 'Phase', 'Group', or 'Ellip'(ticity)


%backwards compatibility
if strcmp(dtype,'PsRF')
    parsed_dtype = {'BW','Ps','def','def'};
    return
elseif strcmp(dtype,'SpRF')
    parsed_dtype = {'BW','Sp','def','def'};
    return
elseif strcmp(dtype,'PsRF_lo')
    parsed_dtype = {'BW','Ps','def','lo'};
    return
elseif strcmp(dtype,'SpRF_lo')
    parsed_dtype = {'BW','Sp','def','lo'};
    return
elseif strcmp(dtype,'SW')
    parsed_dtype = {'SW','Ray','phV',''};
    return
end




dlm = [0,strfind(dtype,'_'),length(dtype)+1]';
pdt = cell(1,length(dlm)-1);
for iparse = 1:length(pdt)
    pdt{iparse} = dtype(dlm(iparse)+1:dlm(iparse+1)-1);
end


if strcmp(pdt{1},'BW')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'BW';
    parsed_dtype{2} = '';
    parsed_dtype{3} = 'def';
    parsed_dtype{4} = 'def';
    % PHASE
    if any(strcmp(pdt,'Ps'))
        parsed_dtype{2} = 'Ps';
    elseif any(strcmp(pdt,'Sp'))
        parsed_dtype{2} = 'Sp'; 
    end
    % TWIN
    if any(strcmp(pdt,'cms'))
        parsed_dtype{3} = 'cms';
    end        
    % FILTF
    if any(strcmp(pdt,'lo'))
        parsed_dtype{4} = 'lo';
    end        

elseif strcmp(pdt{1},'SW')
    % defaults
    parsed_dtype = cell(1,4); 
    parsed_dtype{1} = 'SW';
    parsed_dtype{2} = 'Ray';
    parsed_dtype{3} = 'phV';
    % WAVETYPE
    if any(strcmp(pdt,'Ray'))
        parsed_dtype{2} = 'Ray';
    elseif any(strcmp(pdt,'Lov'))
        parsed_dtype{2} = 'Lov'; 
    elseif any(strcmp(pdt,'HV'))
        parsed_dtype{2} = 'HV'; 
        parsed_dtype{3} = 'HVr';
    end
    % VELTYPE
    if any(strcmp(pdt,'ph'))
        parsed_dtype{3} = 'phV';
    elseif any(strcmp(pdt,'gr'))
        parsed_dtype{3} = 'grV'; 
    end
end






