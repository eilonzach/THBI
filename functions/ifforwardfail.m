function [ iffail ] = ifforwardfail( predata,par )
% [ iffail ] = ifforwardfail( predata, par)
% 
%  quick function to return a true/false on the forward model failing (as
%  described by the predicted data having something horrible in it...)

iffail = false;

if ~isempty(predata.PsRF) 
    if any([predata.PsRF.nsamp]<[predata.PsRF.samprate]*diff(par.datprocess.Twin.PsRF))
        fprintf('Not enough P data!\n');
        iffail = true;
    end
    
    if any(any(isnan([predata.PsRF.ZRT]))) 
        fprintf('NaN P DATA!\n')
        iffail = true;
    end
end

if ~isempty(predata.SpRF)
    if any([predata.SpRF.nsamp]<[predata.SpRF.samprate]*diff(par.datprocess.Twin.SpRF))
        fprintf('Not enough S data!\n');
        iffail = true;
    end

    if any(any(isnan([predata.SpRF.ZRT])))
        fprintf('inhomogeneous!\n'); 
        iffail = true;
    end
end        



end

