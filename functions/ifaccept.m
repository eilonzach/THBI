function accept = ifaccept( E, minE )
%accept = ifaccept( E, minE )
%   function to test whether a new error low enough to be accepted

accept = false;

acc_term = min([1,exp(-0.5*(E-minE))]);
r_deviate = random('unif',0,1,1);

if r_deviate<=acc_term
    accept=true;
end

end

