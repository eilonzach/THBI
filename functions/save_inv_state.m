function save_inv_state(resdir,chainstr,allmodels,misfits) %#ok<INUSD>
% save_inv_state(resdir,chainstr,allmodels,misfits)
%   function to save the state of the inversion into actual .mat files in
%   case the inversion crashes out; we can recover the chains to some
%   extent
    
    save([resdir,'/',chainstr,'_invState'],'allmodels','misfits')

end

