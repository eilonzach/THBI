function [ dat_out ] = clip_main_arrival( dat_in,tt_in,minperiod,PS )
%[ dat_out ] = clip_main_arrival( dat_in,tt_in,minperiod,PS )
%  
% function to clip and taper a data series to get minimise the main arrival
% by locating the centre of the main arrival (the max amplitude point
% across all of the channels) and beginning a taper here, of width half the
% minimum period - this should immediately cut off half the main arrival,
% and taper the remaining half down considerably (leaving about a quarter
% of the main arrival).
% The goal here is to stop the main arrival dominating the cross
% convolution misfit, while not entirely eliminating it or eliminating any
% energy from conversions that come in overlapping/just after the main
% arrival.
% Depending on P/S, we cut from/until this main arrival, respectively.

%% flip direction depending on Ps/Sp
if strcmp(PS,'S')
    dat_in = flipud(dat_in);
    tt_in = flipud(-tt_in);
end

%% find absolute max
[mx,inds]= maxab(dat_in); % max and index
[~,ipp] = maxab(mx(:)); % parent component
amx = mx(ipp);
imx = inds(ipp);
% find opposite swing max
[mx2,imx2] = max(-sign(amx)*dat_in(:,ipp));

%% time of main arrival 
tt_main = tt_in(imx); % default is use main arrival time 

if 2*abs(mx2) > abs(amx) % secondary arrival at least half amp of main arr
%     if strcmp(PS,'P')
        tt_main = tt_in(max([imx,imx2]));% if P, use later time
%     elseif strcmp(PS,'S')
%         tt_main = tt_in(min([imx,imx2]));% if S, use earlier time
%     end
end

%% make one-sided tapered boxcar
win = ones(size(tt_in));
win(tt_in<tt_main) = 0;
% choose the taper indices: [main_arrival to main_arrival + minperiod/2]
itpr = [find(tt_in>=tt_main,1,'first'):find(tt_in<=(tt_main+minperiod/2),1,'last')];
taperlength = length(itpr);
tpr = hanning(2*taperlength);
win(itpr) = tpr(1:taperlength);

%% taper the data
dat_out = dat_in.*(win*ones(1,size(dat_in,2)));


%% reflip
if strcmp(PS,'S')
    dat_out = flipud(dat_out);
end


end

