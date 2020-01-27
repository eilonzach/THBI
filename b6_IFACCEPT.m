function [ ifaccept ] = b6_IFACCEPT( log_likelihood,last_loglikelihood,temp,p_bd,Pm_prior,last_Pm_prior )
% [ ifaccept ] = b6_IFACCEPT( log_likelihood,last_logL,temp,p_bd,Pm_prior,last_Pm_prior  )
% 
% Function to determine whether or not to accept the new model, based on
% its likelihood in comparison to previous likelihoods in the accepted
% model space. 

if nargin<4
    % p_bd is a special term only relevant for birth/death steps
    % includes the probability of the birth/death bit to maintain detailed
    % balance
    p_bd=1;
end

if p_bd == 0
    ifaccept = false;
    return
end

% definitely accept if better (with caveats for prior prob)
if log_likelihood + log(Pm_prior) + log(p_bd) >= last_loglikelihood + log(last_Pm_prior)
    ifaccept = true; 
else % accept according to Metropolis-Hastings law
    r = random('unif',0,1,1);
    ProbRatio = temp * exp(log_likelihood - last_loglikelihood)*(Pm_prior/last_Pm_prior) * p_bd; 
    if ProbRatio >= r
        ifaccept = true;
    else 
        ifaccept = false;
    end

end
    


end

