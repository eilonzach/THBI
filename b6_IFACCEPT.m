function [ ifaccept ] = b6_IFACCEPT( log_likelihood,last_logL,temp,p_bd )
% [ ifaccept ] = b6_IFACCEPT( log_likelihood,last_logL,temp,p_bd )
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

% definitely accept if better
if log_likelihood + log(p_bd) > last_logL
    ifaccept = true; 
else % accept according to Metropolis-Hastings law
    r = random('unif',0,1,1);
    Lr = p_bd*temp*exp(log_likelihood - last_logL); 
    if Lr >= r
        ifaccept = true;
    else 
        ifaccept = false;
    end

end
    


end

