function [log_likelihood_reset,misfit_reset] = b8_LIKELIHOOD_RESET(par,predata,trudata,Kbase,datahparm)
% [Kbase] = make_allkernels(model,Kbase,periods,ID,par)
%   Function to make all surface wave kernels

ifplot = 0;

SWwt = make_SW_weight( par,Kbase,trudata );

misfit_reset = b4_CALC_MISFIT( trudata,predata,par,ifplot,SWwt ); % misfit has structures of summed errors

[ log_likelihood_reset,misfit_reset ] = b5_CALC_LIKELIHOOD( misfit_reset,trudata,datahparm,par);

end

