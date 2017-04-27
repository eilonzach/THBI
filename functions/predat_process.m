function [ data ] = predat_process( data,dtype,par )
%[ dat_out ] = predat_process( dat_in,cp,ifnorm,ifdec )
%   
%  Function to process predicted data by cutting, tapering, and normalising
%  to unit energy. Then decimates to 2*Nyquist.
% 
% INPUT:
%  dat_in  - structure with fields:
%             .ZRT - data in columns, ordered Z,R,T
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

if strcmp(dtype,'SW'), return; end

dat_in = data.(dtype);

dat_out = dat_in;

% loop on number of data time series (each will be 3 component)
for itr = 1:length(dat_in)
if isempty(dat_in(itr).ZRT), continue; end
    
%% make cp
cp = struct('samprate',dat_in(itr).samprate,                 ...
            'pretime',-dat_out(itr).tt(1),  ...
            'prex',-par.datprocess.Twin.(dtype)(1),     ...
            'postx',par.datprocess.Twin.(dtype)(2),     ...                 
            'fhi',par.datprocess.filtf.(dtype)(1),      ...
            'flo',par.datprocess.filtf.(dtype)(2),      ...
            'taperx',0.06,'npoles',2,'norm',0           );


%% clean, filter, window, taper
[ dat_out(itr).ZRT,~,~,~,~,dat_out(itr).tt,~ ] = data_clean(dat_out(itr).ZRT,cp);
% apply windowing (i.e. cut out the zeros from existing taper)
gdind = find( (dat_out(itr).tt>= par.datprocess.Twin.(dtype)(1)) & (dat_out(itr).tt < par.datprocess.Twin.(dtype)(2)) );
dat_out(itr).ZRT = dat_out(itr).ZRT(gdind,:);
dat_out(itr).tt = dat_out(itr).tt(gdind,:);
dat_out(itr).nsamp = length(dat_out(itr).tt);

%% cut some of main arrival
% inwind = (tt_ps >= par.datprocess.Twin.PsRF(1)) & (tt_ps <= par.datprocess.Twin.PsRF(2)); 
% % crop
% predat_ps = predat_ps(inwind,:);
% tt_ps = tt_ps(inwind);

if par.datprocess.clipmain 
    [ dat_out(itr).ZRT ] = clip_main_arrival( dat_out(itr).ZRT,dat_out(itr).tt,1./cp.fhi,dtype(1) );
end

%% normalise to unit energy
if par.datprocess.normdata && ~isempty(dat_out(itr).ZRT)
    normf_ps = dat_out(itr).ZRT(:,1)'*dat_out(itr).ZRT(:,1) + ...
               dat_out(itr).ZRT(:,2)'*dat_out(itr).ZRT(:,2) + ...
               dat_out(itr).ZRT(:,3)'*dat_out(itr).ZRT(:,3);
    dat_out(itr).ZRT = dat_out(itr).ZRT/sqrt(normf_ps);
end
    
%% decimate
if par.datprocess.decdata && ~isempty(dat_out(itr).ZRT)
    resamprate = cp.fhi*4;
    tt_new = [dat_out(itr).tt(1):1./resamprate:dat_out(itr).tt(end)]';
    dat_out(itr).ZRT = interp1(dat_out(itr).tt,dat_out(itr).ZRT,tt_new);
    dat_out(itr).tt = tt_new;
    dat_out(itr).nsamp = length(dat_out(itr).tt);
    dat_out(itr).samprate=resamprate;
end

end % end loop on number of time series

%% put back in
data.(dtype) = dat_out;

end

