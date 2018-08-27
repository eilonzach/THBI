function c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain) %#ok<INUSD>
%  c0_SAVE_OUTPUT(resdir,misfits_perchain,allmodels_perchain)
% 
%  Function to save the important outputs from the inverstion into the
%  resuts directory

save([resdir,'/misfits_perchain_orig'],'misfits_perchain');
save([resdir,'/allmodels_perchain_orig'],'allmodels_perchain');

end
