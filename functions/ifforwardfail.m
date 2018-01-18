function [ iffail ] = ifforwardfail( predata,par )
% [ iffail ] = ifforwardfail( predata, par)
% 
%  quick function to return a true/false on the forward model failing (as
%  described by the predicted data having something horrible in it...)

iffail = false;

if isfield(predata,'BW_Ps') &&~isempty(predata.BW_Ps) 
    if any([predata.BW_Ps.nsamp]<[predata.BW_Ps.samprate]*diff(par.datprocess.Ps.Twin.def))
        fprintf('Not enough P data!\n');
        iffail = true;
    end
    
    if any(any(isnan([predata.BW_Ps.PSV]))) 
        fprintf('NaN P DATA!\n')
        iffail = true;
    end
end

if isfield(predata,'BW_Sp') && ~isempty(predata.BW_Sp)
    if any([predata.BW_Sp.nsamp]<[predata.BW_Sp.samprate]*diff(par.datprocess.Sp.Twin.def))
        fprintf('Not enough S data!\n');
        iffail = true;
    end

    if any(any(isnan([predata.BW_Sp.PSV])))
        fprintf('inhomogeneous!\n'); 
        iffail = true;
    end
end        

end

