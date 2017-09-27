function Vp  = mantle_vs2vp( Vs,Zkm )
% Vp  = mantle_vs2vp( Vs,Zkm )
%   Scaling of Vs to Vp for mantle rocks, using Vp/Vs ratio at each depth
%   from AK135.


% akmod = ak135('depths',Zkm,'crust',false);

try
%     Vp = Vs.*(akmod.vp./akmod.vs);
    Vp = Vs.*ak135vpvs(Zkm);
catch
    Vp=1.81*Vs
    error
end

end

function vpvs = ak135vpvs(Zkm)

% values are taken from ak135 function with no crust
% columns are Z, vpvs_nocrust
a = [
     0.000      1.7970
    20.000      1.7957
    20.0001     1.7957
    35.000      1.7946
    35.0001     1.7946
    77.500      1.7918
   120.000      1.7889
   165.000      1.8130
   210.000      1.8371
   210.0001     1.8351
   260.000      1.8404
   310.000      1.8452
   360.000      1.8498
   410.000      1.8425];

vpvs = interp1(a(:,1),a(:,2),Zkm);

end
