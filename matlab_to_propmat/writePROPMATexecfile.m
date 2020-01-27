function writePROPMATexecfile( execfile, modfile,ifile,ofile0,ofile1,ofile2,oimfile,odatfile,inc,PS,soT,obsdist,ocomps)
%writePROPMATexecfile( execfile, modfile,ifile,ofile0,ofile1,ofile2,oimfile,odatfile,inc,PS,soT,obsdist,ocomps)
%   
% Function to write execution file for PropMatrix code
% 
% INPUTS:
%  execfile  - name of execution file to write
%  modfile   - file with model description
%  ifile     - name of synth.in file
%  ofile0    - name of synth.out file, not used again
%  ofile1    - name of 1st file with synth output, for sourc1 to use
%  ofile2    - name of 2nd file with synth output, for sourc1 to use
%  oimfile   - name of imaginary output file, will be sour1 output
%  odatfile  - name of file with output data from sourc1
%  inc       - incidence angle
%  PS        - choice of Ps or Sp (strings)
%  soT       - source wavelet period (s)
%  obsdist   - observation distance (km)
%  ocomps    - output format: 1 = x,y,z  or 2 = R,T,U 


if nargin < 10
    PS = 'Ps';
end
if nargin < 11
    soT = 1;
end
if nargin < 12 
    obsdist = 0;
end
if nargin < 13
    ocomps = 2;
end
    
if strcmp(PS,'Ps')
    excite = '1 0 0';
elseif strcmp(PS,'Sp')
    excite = '0 1 0';
end

if exist(execfile,'file')==2
    delete(execfile); % kill if it is there 
end

%% write synth.in parameter file
fid = fopen(execfile,'w');
fprintf(fid,'#!/bin/csh\n');
fprintf(fid,'set xdir=/Users/zeilon/Work/codes/PropMatrix/bin\n');
fprintf(fid,'$xdir/synth <<!\n');
fprintf(fid,'%s\n',modfile);
fprintf(fid,'%s\n',ofile0);
fprintf(fid,'%s\n',ofile1);
fprintf(fid,'%s\n',ofile2);
fprintf(fid,'%s\n',ifile);
fprintf(fid,'%.1f\n',inc);
fprintf(fid,'!\n');
%
fprintf(fid,'echo "Done propagation, now convolving source"\n');
%
fprintf(fid,'$xdir/sourc1 <<!\n');
fprintf(fid,'%s\n',odatfile);
fprintf(fid,'%s\n',ofile1);
fprintf(fid,'%s\n',ofile2);
fprintf(fid,'%s\n',oimfile);
fprintf(fid,'%s\n',excite);
fprintf(fid,'%.1f\n',soT);
fprintf(fid,'%.1f\n',obsdist);
fprintf(fid,'%.0f\n',ocomps);
fprintf(fid,'!\n');
%
fprintf(fid,'echo "Done source convolution, cleaning up..."\n');
%
% fprintf(fid,'rm synth.out3*\n');
fprintf(fid,'rm %s\n',oimfile);

fclose(fid);

end




