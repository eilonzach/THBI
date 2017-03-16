function ikernelfiles = writeKERNELCALCexecfile(swperiods,execfile,stripfile,eigfile,qmod,tabfile,qfile,kernelfile,ikprefix,logfile)
% ikernelfiles = writeKERNELCALCexecfile(swperiods,execfile,stripfile,eigfile,tabfile,qfile,kernelfile,ikprefix,logfile)
%   
% Function to write execution file for kernel calculator having run MINEOS code
% 
% INPUTS:
%  swperiods - vector of surface wave periods (will be rounded to integers)
%  execfile  - name of executable file 
%  modefile  - name of strip file
%  eigfile   - name of eigenfunctions output binary file
%  qmod      - name of qmod file with details about (unused) Q model
%  tabfile   - name of output table file
%  qfile     - name of output q file (= output from mineos_q; not the same as the Q model)
%  kernelfile- name of output frechet kernel file -- this is a key file
%  ikprefix  - prefix for output individual kernel files at each freq.
%  logfile   - name of file to print screen output to

ikernelfiles = cell({});

if exist(execfile,'file')==2
    delete(execfile); % kill if it is there 
end

%% write synth.in parameter file
fid = fopen(execfile,'w');
fprintf(fid,'#!/bin/csh\n');
%
fprintf(fid,'#\n');
fprintf(fid,'set xdir=/Users/zeilon/Work/codes/MINEOS_Dalton/bin\n');
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Stripping mineos"\n');
%
fprintf(fid,'#\n');
%
fprintf(fid,'$xdir/mineos_strip <<! > %s\n',logfile);
fprintf(fid,'%s\n',stripfile);
fprintf(fid,'%s\n',eigfile);
fprintf(fid,'\n');
fprintf(fid,'!\n');
%
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Done stripping, now calculating tables"\n');
%
fprintf(fid,'#\n');
%
fprintf(fid,'$xdir/mineos_table <<! >> %s\n',logfile);
fprintf(fid,'%s\n',tabfile);
fprintf(fid,'40000\n');
fprintf(fid,'0 200.1\n');
fprintf(fid,'0 3000\n');
fprintf(fid,'%s\n',qfile);
fprintf(fid,'%s\n',stripfile);
fprintf(fid,'\n');
fprintf(fid,'!\n');
%
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Creating branch file"\n');
%
fprintf(fid,'#\n');
fprintf(fid,'# to create branch file needed for frechet derivatives:\n');
fprintf(fid,'# second line says stop searching (or could add more parameters to search)\n');
fprintf(fid,'# 3rd line gives frequency range to search in (mHz)\n');
fprintf(fid,'#\n');
%
fprintf(fid,'$xdir/plot_wk <<! > %s\n',logfile);
fprintf(fid,'table %s_hdr\n',tabfile);
fprintf(fid,'search\n');
fprintf(fid,'1 0.0 200.05\n');
fprintf(fid,'99 0 0\n');
fprintf(fid,'branch\n');
fprintf(fid,'\n');
fprintf(fid,'quit\n');
fprintf(fid,'!\n');
%
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Making frechet kernels binary"\n');
%
fprintf(fid,'#\n');
%
fprintf(fid,'trash %s\n',kernelfile);
fprintf(fid,'$xdir/frechet_cv <<! >> %s\n',logfile);
fprintf(fid,'%s\n',qmod);
fprintf(fid,'%s_hdr.branch\n',tabfile);
fprintf(fid,'%s\n',kernelfile);
fprintf(fid,'%s\n',eigfile);
fprintf(fid,'0\n');
fprintf(fid,'\n');
fprintf(fid,'!\n');
%
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Writing kernel files for each period"\n');
for ip = 1:length(swperiods)
%
fprintf(fid,'#\n');
%
fprintf(fid,'$xdir/draw_frechet_gv <<!\n');
fprintf(fid,'%s\n',kernelfile);
fprintf(fid,'%s_cvfrechet_%.0fs\n',ikprefix,round(swperiods(ip)));
fprintf(fid,'%.0f\n',round(swperiods(ip)));
fprintf(fid,'!\n');
%
ikernelfiles{ip} = sprintf('%s_cvfrechet_%.0fs',ikprefix,round(swperiods(ip)));
end
fprintf(fid,'#\n');
%% =======================================================================
fprintf(fid,'echo "Done velocity calculation, cleaning up..."\n');
%
% fprintf(fid,'rm synth.out3*\n');
fprintf(fid,'rm %s\n',logfile);

fclose(fid);

end




