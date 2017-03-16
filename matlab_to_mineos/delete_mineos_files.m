function delete_mineos_files( ID,swperiods )
%delete_mineos_files( ID )
%   Function to delete all files associted with mineos running

%% filenames
if ~ischar(ID), ID = num2str(ID);end
execfile = [ID,'.run_mineos'];
cardfile = [ID,'.model'];
eigfile = [ID,'.eig'];
ofile1 = [ID,'.asc1'];
logfile = [ID,'.log'];
execfile_k = [ID,'.run_kernels'];
stripfile = [ID,'.strip'];
tabfile = [ID,'.table'];
qfile = [ID,'.q'];
kernelfile = [ID,'.cvfrechet'];

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

