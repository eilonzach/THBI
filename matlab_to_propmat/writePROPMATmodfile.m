function writePROPMATmodfile( model1D,modfile )
%writePROPMATmodfile( model1D,modfile )
%   
% Function to write model file for PropMatrix code
% 
% INPUTS:
%  model1D  - model to be written, with fields comprising [nlay*1] vectors:
%   .zlayt  - depths to layer tops (m)
%   .zlayb  - depths to layer bottoms (m)
%   .Vp     - Vp (m/s)
%   .Vs     - Vs (m/s)
%   .rho    - density (kg/m^3)
%   .xi     - S radial anisotropy (Vsh/Vsv)^2
%   .phi    - S radial anisotropy (Vpv/Vph)^2
%   .eta    - eta factor, F/(A-2L)
%   .VpAnis - Vp azimuthal anisotropy percent (peak-to-peak)
%   .VsAnis - Vs azimuthal anisotropy percent (peak-to-peak)
%   .faz    - azimuth of fast axis of anisotropy (deg from N)
%   .fpl    - plunge of fast axis of anisotropy (deg from horiz)
%   .intstr - strike of dipping interface (deg from N)
%   .intdip - dip of dipping interface (deg from horiz)
%  modfile  - file name for output file


nlay = length(model1D.zlayt);
% defaults
wmod = struct('thick',model1D.zlayb(:)-model1D.zlayt(:),...
              'rho',nan(nlay,1),'alph',nan(nlay,1),'beta',nan(nlay,1),...
              'xi',ones(nlay,1),'phi',ones(nlay,1),'eta',ones(nlay,1),...
              'iso',ones(nlay,1),'Panis',zeros(nlay,1),'Sanis',zeros(nlay,1),...
              'tr',ones(nlay,1),'pl',zeros(nlay,1),...
              'strk',zeros(nlay,1),'dip',zeros(nlay,1));          
% infinite half space at bottom
          
modparms = fieldnames(model1D);
for ii = 1:length(modparms)
    switch modparms{ii}
        case {'Vp'}
            wmod.alph = model1D.(modparms{ii});
        case {'Vs'}
            wmod.beta = model1D.(modparms{ii});
        case {'rho','density'}
            wmod.rho = model1D.(modparms{ii});
        case {'Panis','Vp_anis','anisP'}
            wmod.Panis = model1D.(modparms{ii});
        case {'Sanis','Vs_anis','anisS'}
            wmod.Sanis = model1D.(modparms{ii});    
        case {'xi','Xi','VshVsv'}
            wmod.xi = model1D.(modparms{ii});    
        case {'phi','Phi','VpvVph'}
            wmod.phi = model1D.(modparms{ii});    
        case {'eta','Eta'}
            wmod.eta = model1D.(modparms{ii});    
    end
end


%% Radial anisotropy
[ ys,xs ] = VsvVsh_from_VsXi( wmod.beta,model1D.xi );
[ yp,xp ] = VpvVph_from_VpPhi( wmod.alph,model1D.phi );

anisflag = zeros(nlay,1);
for kk = 1:nlay
    if all((xs==ys) & (xp==yp) & (wmod.eta==ones(nlay,1))) 
        anisflag(kk) = 1;
    else
        anisflag(kk) = 1;
    end
end
A = wmod.rho .* xp.^2; % if no anis, xp just P-wave modulus
C = wmod.rho .* yp.^2; % if no anis, yp just P-wave modulus
L = wmod.rho .* ys.^2; % if no anis, ys just shear modulus
N = wmod.rho .* xs.^2; % if no anis, xs just shear modulus
F = wmod.eta.*(A - 2*L);     % if no anis, is just lambda

%% write model  file
fid = fopen(modfile,'w');
for kk = 1:nlay
%     fprintf(fid,'%6.0f %6.1f %6.1f %6.1f  %1.0f  %4.2f %4.2f %3.0f %2.0f %3.0f %2.0f\n',...
%         1e3*wmod.thick(kk),1e3*wmod.rho(kk),1e3*wmod.alph(kk),1e3*wmod.beta(kk),wmod.iso(kk),wmod.Panis(kk),wmod.Sanis(kk),wmod.tr(kk),wmod.pl(kk),wmod.strk(kk),wmod.dip(kk));
    
    % data line
    fprintf(fid,'layer %.0f\n',kk);                                    % layer description
    fprintf(fid,'%-2.0f %4.2f %6.1f\n',anisflag(kk),wmod.rho(kk),wmod.thick(kk)); % layer# density thickness(km)
    fprintf(fid,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',A(kk),A(kk)-2*N(kk), F(kk),0,0,0);
    fprintf(fid,'        %7.3f %7.3f %7.3f %7.3f %7.3f\n',A(kk),F(kk),0,0,0);
    fprintf(fid,'                %7.3f %7.3f %7.3f %7.3f\n',C(kk),0,0,0);
    fprintf(fid,'                        %7.3f %7.3f %7.3f\n',L(kk),0,0);
    fprintf(fid,'                                %7.3f %7.3f\n',L(kk),0);
    fprintf(fid,'                                        %7.3f\n',N(kk));
end
% half space at the bottom
fprintf(fid,'layer halfspace\n');                                    % layer description
fprintf(fid,'%-2.0f %4.2f %6.1f\n',anisflag(nlay),wmod.rho(nlay),1); % layer# density thickness(km)
fprintf(fid,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n',A(nlay),A(nlay)-2*N(nlay), F(nlay),0,0,0);
fprintf(fid,'        %7.3f %7.3f %7.3f %7.3f %7.3f\n',A(nlay),F(nlay),0,0,0);
fprintf(fid,'                %7.3f %7.3f %7.3f %7.3f\n',C(nlay),0,0,0);
fprintf(fid,'                        %7.3f %7.3f %7.3f\n',L(nlay),0,0);
fprintf(fid,'                                %7.3f %7.3f\n',L(nlay),0);
fprintf(fid,'                                        %7.3f\n',N(nlay));
    
    
fclose(fid);


end




