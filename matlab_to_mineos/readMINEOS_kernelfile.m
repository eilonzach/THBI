function kernel = readMINEOS_kernelfile(kernelfile)
% kernel = readMINEOS_kernelfile(kernelfile)
%  
%  Function to read MINEOS kernel file to get the perturbational
%  sensitivity kernels for dc/c as a function of dvpv/vpv, dvsv/vsv etc

fid = fopen(kernelfile,'r');
C = textscan(fid,'%f %f %f %f %f %f %f'); % Z,Vsv, Vpv, Vsh, Vph, eta, rho
fclose(fid);

kernel = struct('Z',6371e3-flipud(C{1}),... % flipuds so Z=0 at the top
                'Vsv',flipud(C{2}),'Vsh',flipud(C{4}),...
                'Vpv',flipud(C{3}),'Vph',flipud(C{5}),...
                'eta',flipud(C{6}),'rho',flipud(C{7}));
end

