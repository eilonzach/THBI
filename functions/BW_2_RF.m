function [trudata_rf] = bw2rf(trudata,par)
%trudata_rf] = bw2rf(trudata,RFparms)
%   Function to take body wave data and swap the traces for receiver
%   functions, where the parent phase is deconvolved from a) the parent
%   (i.e. autodeconvolved) to and b) the daughter, to make new parent and
%   daughter traces

RFparms = par.synth.RFparms;
trudata_rf = trudata;

% Find the body wave data types to make RFs for
trudtypes = fieldnames(trudata);
for idt = 1:length(trudtypes)    

    dtype = trudtypes{idt};
	RFdtype = regexprep(dtype,'BW','RF');
 
    pdt = parse_dtype( trudtypes{idt} );
    if strcmp(pdt{1},'BW') && any(ismember(par.inv.datatypes,RFdtype))
        % Identify parent channel (assume that P means parent is 1, S means
        % parent is 2!)
        ip = find(strcmp(pdt{2}(1),{'P','S'})); 
        % daughter channel is whichever of those is not the parent...
        id = setdiff([1;2],ip);
        
        % grab data
        Parent = trudata.(trudtypes{idt}).PSV(:,ip);
        Daughter = trudata.(trudtypes{idt}).PSV(:,id);
      
        
        %---RF parameters---%
    %=================================%
    %input: time domain deconvolution
%     RFparms.gauss_t = 0.5;
%     RFparms.accept_mis=1e-10; %accepted misfit
%     RFparms.itmax=50; %number of iterations
% 
%     %input: freq domain deconvolution
%     RFparms.TB = 1.5;%1.5; %period*bandwidth
%     RFparms.NT = 2;%2; %number of tapers
%     RFparms.tag = 'synth'; %data or synth (synthetic) - Zach uses 'synth' & this gives better results
%     RFparms.Poverlap = 0.95; %fractional window overlap
    %=================================%
        
        % Autodeconvolve parent
        if strcmp(RFparms.method,'ETMTM')
            
            [tt_RFp, RFp] = ETMTM(Parent',Parent',pdt{2},... % Parent, Daughter, 'Ps'/'Sp'
                                RFparms.TB,RFparms.NT,RFparms.tag,... % period for gaussian, number of tapers, 'synth'/'data' (synth seems to work better?)
                                1./trudata.(trudtypes{idt}).samprate,... % dt
                                trudata.(trudtypes{idt}).nsamp./trudata.(trudtypes{idt}).samprate,... % wlength
                                RFparms.Poverlap);
        % Deconvolve daughter
            [tt_RFd, RFd] = ETMTM(Parent',Daughter',pdt{2},... % Parent, Daughter, 'Ps'/'Sp'
                                RFparms.TB,RFparms.NT,RFparms.tag,... % period for gaussian, number of tapers, 'synth'/'data' (synth seems to work better?)
                                1./trudata.(trudtypes{idt}).samprate,... % dt
                                trudata.(trudtypes{idt}).nsamp./trudata.(trudtypes{idt}).samprate,... % wlength
                                RFparms.Poverlap);
        
        elseif strcmp(RFparms.method,'IDRF')
            
            [tt_RFp, RFp]=IDRF(Parent',Parent',...
                           1./trudata.(trudtypes{idt}).samprate,... % dt
                           minmax(trudata.(trudtypes{idt}).tt),RFparms.gauss_t,...
                           RFparms.accept_mis,RFparms.itmax);
            [tt_RFd, RFd]=IDRF(Parent',Daughter',...
                           1./trudata.(trudtypes{idt}).samprate,... % dt
                           minmax(trudata.(trudtypes{idt}).tt),RFparms.gauss_t,...
                           RFparms.accept_mis,RFparms.itmax);
                            
        end
                            
%         figure(21); clf; subplot(211);plot(tt_RFp,RFp);subplot(212);plot(tt_RFd,RFd)
        
        %% insert back into trudata - interpolating times onto original time
        trudata_rf.(trudtypes{idt}).PSV(:,ip) = interp1(tt_RFp,RFp,trudata.(trudtypes{idt}).tt).*maxab(maxab(trudata.(trudtypes{idt}).PSV));
        trudata_rf.(trudtypes{idt}).PSV(:,id) = interp1(tt_RFd,RFd,trudata.(trudtypes{idt}).tt).*maxab(maxab(trudata.(trudtypes{idt}).PSV));

        % rename data type to RF
        trudata_rf.(RFdtype) = trudata_rf.(dtype);
        if ~any(ismember(par.inv.datatypes,dtype))
            trudata_rf = rmfield(trudata_rf,dtype);
        end

%         figure; 
%         subplot(211);
%         plot(trudata.(trudtypes{idt}).tt,[trudata.(trudtypes{idt}).PSV(:,ip),trudata.(trudtypes{idt}).RF_PSV(:,ip)])
%         subplot(212);
%         plot(trudata.(trudtypes{idt}).tt,[trudata.(trudtypes{idt}).PSV(:,id),trudata.(trudtypes{idt}).RF_PSV(:,id)])
        
        
    end
end


end

