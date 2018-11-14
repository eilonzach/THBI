function [zlayt,zlayb,vslay,z,vs,vp,rho] = RD_1D_Vprofile
% [zlayt,zlayb,vslay,z,vs] = RD_1D_Vprofile
%   Function to grab the calculated 1D shear velocity profile from a
%   digitised version of the blue curve in Riddhi Dave and Aibing Li's 2016
%   paper's online supplemental figure DR3_B

a = load('~/Dropbox/Dave_Li_phV/Wyoming1D.txt','-ascii');

z = a(:,2);
vs = a(:,1);

zb = z(2:2:end);
zt = z(1:2:end);
zdisc = round_level(mean([zb(1:end-1),zt(2:end)],2),5);
zlayb = [zdisc;300];
zlayt = [0;zdisc];

vslay = round_level(mean([vs(1:2:end),vs(2:2:end)],2),0.005);

nlay = length(vslay);

z = reshape([zlayt';zlayb'],2*nlay,1);
vs = reshape([vslay';vslay'],2*nlay,1);
% figure, plot(vs,z);set(gca,'ydir','reverse','xlim',[3.2 4.8])

zmoh = z(find(vs>4.1,1,'first'));
vp = 1.81* vs;
rho(z<zmoh) = sed_vs2rho(vs(z<zmoh));
rho(z>=zmoh) = mantle_vs2rho(vs(z>=zmoh),z(z>=zmoh));
rho = rho(:);


end

