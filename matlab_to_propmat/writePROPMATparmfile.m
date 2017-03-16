function writePROPMATparmfile( parmfile, Vinc, Nlay,nsamps,samprate,cutf)
%writePROPMATparmfile( parmfile, Vinc, Nlay,nsamps,samprate)
%   
% Function to write synth.in file for PropMatrix code
% 
% INPUTS:
%  parmfile  - name of output parameter file to write
%  Vinc      - velocity of incident wave - should be same as bottom layer
%              of model
%  Nlay      - number of layers in model
%  nsamps    - number of samples in output time series (will round up to
%              nearest power of 2)
%  samprate  - sample rate of output time series
%  cutf      = cutting freuquency

%% write synth.in parameter file
fid = fopen(parmfile,'w');
fprintf(fid,'%-10.2f(velocity  of incident wave)\n',Vinc);
fprintf(fid,'%-10.0f(no. of layers, including half-space)\n',Nlay);
fprintf(fid,'%-10.0f(no. of FFT)\n',nsamps);
fprintf(fid,'%-10.0f(no. of points per second, npps)\n',samprate);
fprintf(fid,'%.0f       (cutting frequency)\n',cutf);
fprintf(fid,'1         (2=bottom motion, 1= surface motion)\n');
fclose(fid);

end




