function writeRAYSUMmodfile( model1D,modfile )
%writeRAYSUMmodfile( model1D,modfile )
%   
% Function to write model file for RAYSUM code
% 
% INPUTS:
%  model1D  - model to be written, with fields comprising [nlay*1] vectors:
%   .zlayt  - depths to layer tops (m)
%   .zlayb  - depths to layer bottoms (m)
%   .Vp     - Vp (m/s)
%   .Vs     - Vs (m/s)
%   .rho    - density (kg/m^3)
%   .VpAnis - Vp anisotropy percent (peak-to-peak)
%   .VsAnis - Vs anisotropy percent (peak-to-peak)
%   .faz    - azimuth of fast axis of anisotropy (deg from N)
%   .fpl    - plunge of fast axis of anisotropy (deg from horiz)
%   .intstr - strike of dipping interface (deg from N)
%   .intdip - dip of dipping interface (deg from horiz)
%  modfile  - file name for output file

nlay = length(model1D.zlayt);
% defaults
wmod = struct('thick',model1D.zlayb(:)-model1D.zlayt(:),...
              'rho',nan(nlay,1),'alph',nan(nlay,1),'beta',nan(nlay,1),...
              'iso',ones(nlay,1),'Panis',zeros(nlay,1),'Sanis',zeros(nlay,1),...
              'tr',ones(nlay,1),'pl',zeros(nlay,1),...
              'strk',zeros(nlay,1),'dip',zeros(nlay,1));          
% infinite half space at bottom
wmod.thick(end)=0;
          
modparms = fieldnames(model1D);
for ii = 1:length(modparms)
    switch modparms{ii}
        case {'Vp','Vpv','vpv'}
            wmod.alph = model1D.(modparms{ii});
        case {'Vs','Vsv','vsv'}
            wmod.beta = model1D.(modparms{ii});
        case {'rho','density'}
            wmod.rho = model1D.(modparms{ii});
        case {'Panis','Vp_anis','anisP'}
            wmod.Panis = model1D.(modparms{ii});
        case {'Sanis','Vs_anis','anisS'}
            wmod.Sanis = model1D.(modparms{ii});    
    end
end





%% write
fid = fopen(modfile,'w');
fprintf(fid,'#thick    rho   alph   beta iso  %%P  %%S    tr pl  st di\n');
for kk = 1:nlay
    fprintf(fid,'%6.0f %6.1f %6.1f %6.1f  %1.0f  %4.2f %4.2f %3.0f %2.0f %3.0f %2.0f\n',...
        1e3*wmod.thick(kk),1e3*wmod.rho(kk),1e3*wmod.alph(kk),1e3*wmod.beta(kk),wmod.iso(kk),wmod.Panis(kk),wmod.Sanis(kk),wmod.tr(kk),wmod.pl(kk),wmod.strk(kk),wmod.dip(kk));
end
fclose(fid);



end

