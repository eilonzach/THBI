function Vp  = mantle_vs2vp( Vs,Zkm )
% Vp  = mantle_vs2vp( Vs,Zkm )
%   Scaling of Vs to Vp for mantle rocks, using Vp/Vs ratio at each depth
%   from AK135.

akmod = ak135('depths',Zkm);

try
    Vp = Vs.*(akmod.vp./akmod.vs);
catch
    Vp=1.81*Vs
    error
end

end

