function [ ifaccept ] = b6_IFACCEPT( log_likelihood,misfits,temp )
% [ ifaccept ] = b6_IFACCEPT( log_likelihood,misfits,temp)
% 
% Function to determine whether or not to accept the new model, based on
% its likelihood in comparison to previous likelihoods in the accepted
% model space. 

if nargin<4
    p_bd=1;
end

% definitely accept if better
if log_likelihood > misfits.lastlogL
    ifaccept = true; 
else % accept according to Metropolis-Hastings law
    r = random('unif',0,1,1);
    Lr = temp*exp(log_likelihood - misfits.lastlogL);
    if Lr >= r
        ifaccept = true;
    else 
        ifaccept = false;
    end

end
    


end

