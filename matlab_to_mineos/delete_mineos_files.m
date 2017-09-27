function delete_mineos_files( ID,R_or_L )
%delete_mineos_files( ID )
%   Function to delete all files associted with mineos running

%% filenames
if ~ischar(ID), ID = num2str(ID);end
ID = [ID,R_or_L];
execfile = [ID,'.run_mineos'];  if exist(execfile,'file')~=2, execfile = []; end
cardfile = [ID,'.model'];       if exist(cardfile,'file')~=2, cardfile = []; end
eigfile = [ID,'.eig'];          if exist(eigfile,'file')~=2, eigfile = []; end
ofile1 = [ID,'.asc1'];          if exist(ofile1,'file')~=2, ofile1 = []; end
logfile = [ID,'.log'];          if exist(logfile,'file')~=2, logfile = []; end
execfile_k = [ID,'.run_kernels'];if exist(execfile_k,'file')~=2, execfile_k = []; end
stripfile = [ID,'.strip'];      if exist(stripfile,'file')~=2, stripfile = []; end
tabfile = [ID,'.table'];        if exist(tabfile,'file')~=2, tabfile = []; end
qfile = [ID,'.q'];              if exist(qfile,'file')~=2, qfile = []; end
kernelfile = [ID,'.cvfrechet']; if exist(kernelfile,'file')~=2, kernelfile = []; end

% preamble
wd = pwd;
cd('/Users/zeilon/Documents/MATLAB/BayesianJointInv/matlab_to_mineos');
%% do the deleting

delete(execfile,cardfile,eigfile,ofile1,qfile);
if exist(execfile_k,'file')==2
    delete(execfile_k,stripfile,tabfile,kernelfile,[tabfile,'_hdr'],[tabfile,'_hdr.branch']);
    for ip = 1:length(swperiods), 
        delete(ikernelfiles{ip}); 
    end
end
if exist(logfile,'file')==2, delete(logfile); end

% postamble
cd(wd);

end

