function [ code ] = dtype_code( dtype )
%[ code ] = dtype_code( dtype )
%   simple function to deliver shorthand for the data type

if strcmp(dtype,'SW')
    code = 'SW';
elseif strcmp(dtype,'PsRF')
    code = 'ps';
elseif strcmp(dtype,'SpRF')
    code = 'sp';
elseif strcmp(dtype,'PsRF_lo')
    code = 'ps_lo';
elseif strcmp(dtype,'SpRF_lo')
    code = 'sp_lo';
end

end

