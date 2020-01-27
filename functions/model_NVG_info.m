function [nvg_z,nvg_w,nvg_a] = model_NVG_info(model)
% [nvg_z,nvg_w,nvg_a] = model_NVG_info(model)
%   Function to get low-velocity zone data for seismic model of the
%   final_model format. The outputs are the depth of the centre (z), width
%   (w) and percentage fractional amplitude (a) of the first nvegative
%   velocity zone, where the amplitude is defined as the fractional
%   difference between the max and min velocities at the top and bottom of
%   the zone, divided by their mean.

%% Find max mantle velocity immediately below moho

z = model.Z;
Vs = model.VSav;
zmoh = model.Zd(2).mu;

% gradient of mantle velocity change
gVs = gradient(Vs,z);

%  point of maximum velocity in 'lithosphere'
maxVz = z(find(gVs<0.00022 & z>(zmoh+5),1,'first')); % lith at least 5 km...
minVz = z(find(gVs>-0.00002 & z>(maxVz+10),1,'first')); % no lith smaller than 10 km...

if isempty(minVz)
    nvg_z = nan;
    nvg_w = nan;
    nvg_a = nan;
    return
end
    

Vs_max = Vs(z==maxVz);
Vs_min = Vs(z==minVz);

nvg_z = mean([maxVz,minVz]);
nvg_w = diff([maxVz,minVz]);
nvg_a = 100*(Vs_max-Vs_min)/mean([Vs_max,Vs_min]);

% figure(45);  clf
% subplot(121),plot(Vs,z,'linewidth',2); set(gca,'ydir','reverse','ylim',[00 200],'xlim',[4.0 4.8])
% subplot(122),plot(gVs,z,'linewidth',2); set(gca,'ydir','reverse','ylim',[00 200],'xlim',[-0.03 0.03])


end

