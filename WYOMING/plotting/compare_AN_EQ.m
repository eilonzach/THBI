% RSSD
slat = 44.121+2;
slon = -104.04-2;

%% Earthquake data:
[ freqs, Eperiods ] = get_freqs;
EphV_period = disp_curve_EQ_latlon(freqs,slat,slon);

%% AN data:
ANperiods = [8:2:32,36,40]';
ANphV_period = disp_curve_AN_latlon(ANperiods,slat,slon);

%% composite
[ periods, phVs ] = disp_curve_latlon(slat,slon,max(ANperiods));

%% plotting
figure(54), clf
plot(periods,phVs,'o-k','linewidth',2)
hold on
plot(Eperiods,EphV_period,'o-b',ANperiods,ANphV_period,'o-r')
