function [ data ] = predat_process( data,dtype,par )
%[ dat_out ] = predat_process( dat_in,cp,ifnorm,ifdec )
%   
%  Function to process predicted data by cutting, tapering, and normalising
%  to unit energy. Then decimates to 2*Nyquist.
% 
% INPUT:
%  dat_in  - structure with fields:
%             .PSV - data in columns, ordered Z,R,T
%             .tt  - tt vector, w/ significant padding before phase arrival
%             .nsamp - number of samples
%  cp      - cleaning_parm structure, as for use with data_clean
%  ifnorm  - option [true/false] to normalise data to unit energy
%  ifdec   - option [true/false] to decimate data to 2*fNyq
%       
% OUTPUT:
%  dat_out - similar structure, but with the data having been cut,
%               tapered with a 5 second hanning window at either end, and
%               normalised to unit energy.
% 
%  Z. Eilon 09/2016

[ pdt ] = parse_dtype( dtype );

if ~strcmp(pdt{1},'BW'), return; end

dat_in = data.(dtype);

dat_out = dat_in;

% loop on number of data time series (each will be 3 component)
for itr = 1:length(dat_in)
if isempty(dat_in(itr).PSV), continue; end

%% cut mis-rotated main arrival out of daughter
if isfield(par.datprocess,'clipdaughter') && par.datprocess.clipdaughter==true; 
    % set daughter channel and clips
    if strcmp(pdt{2}(1),'P')
        chd = 2; winkeep = [0  max(dat_out(itr).tt)+2]; tapertime=1;
    elseif strcmp(pdt{2}(1),'S'), 
        chd = 1; winkeep = [min(dat_out(itr).tt)-2 -0.5]; tapertime=1.5;
    end 
    [ dat_out(itr).PSV(:,chd) ] = flat_hanning_win(dat_out(itr).tt,dat_out(itr).PSV(:,chd),winkeep(1),winkeep(2),tapertime);
end
    
%% make cp
cp = struct('samprate',dat_in(itr).samprate,                 ...
            'pretime',-dat_out(itr).tt(1),  ...
            'prex',-par.datprocess.(pdt{2}).Twin.(pdt{3})(1),     ...
            'postx',par.datprocess.(pdt{2}).Twin.(pdt{3})(2),     ...                 
            'fhi',par.datprocess.(pdt{2}).filtf.(pdt{4})(1),      ...
            'flo',par.datprocess.(pdt{2}).filtf.(pdt{4})(2),      ...
            'taperx',0.1,'npoles',2,'norm',0           );


%% clean, filter, window, taper
[ dat_out(itr).PSV,~,~,~,~,dat_out(itr).tt,~ ] = data_clean(dat_out(itr).PSV,cp);
% apply windowing (i.e. cut out the zeros from existing taper)
gdind = find( (dat_out(itr).tt>= par.datprocess.(pdt{2}).Twin.(pdt{3})(1)) & (dat_out(itr).tt < par.datprocess.(pdt{2}).Twin.(pdt{3})(2)) );
dat_out(itr).PSV = dat_out(itr).PSV(gdind,:);
dat_out(itr).tt = dat_out(itr).tt(gdind,:);
dat_out(itr).nsamp = length(dat_out(itr).tt);

%% cut some of main arrival
if par.datprocess.clipmain 
    [ dat_out(itr).PSV ] = clip_main_arrival( dat_out(itr).PSV,dat_out(itr).tt,1./cp.fhi,pdt{2}(1) );
end
% inwind = (tt_ps >= par.datprocess.Twin.PsRF(1)) & (tt_ps <= par.datprocess.Twin.PsRF(2)); 
% % crop
% predat_ps = predat_ps(inwind,:);
% tt_ps = tt_ps(inwind);


% inwind = (tt_ps >= par.datprocess.Twin.PsRF(1)) & (tt_ps <= par.datprocess.Twin.PsRF(2)); 
% % crop
% predat_ps = predat_ps(inwind,:);
% tt_ps = tt_ps(inwind);


%% normalise to unit energy, flip so max is positive
if par.datprocess.normdata && ~isempty(dat_out(itr).PSV)
    norm2f_ps = dat_out(itr).PSV(:,1)'*dat_out(itr).PSV(:,1) + ...
               dat_out(itr).PSV(:,2)'*dat_out(itr).PSV(:,2);
    dat_out(itr).PSV = dat_out(itr).PSV/sqrt(norm2f_ps)./sign(maxab(dat_out(itr).PSV(:)));
end
    
%% decimate
if par.datprocess.decdata && ~isempty(dat_out(itr).PSV)
    resamprate = cp.fhi*4;
    tt_new = [dat_out(itr).tt(1):1./resamprate:dat_out(itr).tt(end)]';
    dat_out(itr).PSV = interp1(dat_out(itr).tt,dat_out(itr).PSV,tt_new);
    dat_out(itr).tt = tt_new;
    dat_out(itr).nsamp = length(dat_out(itr).tt);
    dat_out(itr).samprate=resamprate;
end

end % end loop on number of time series

%% put back in
data.(dtype) = dat_out;

end

