function [traces,tt,align] = readraysumtr(trfile)
% [traces,tt] = readraysumtr(trfile)

fid=fopen(trfile,'r');

% Read header:
line=fgetl(fid);
while line(1) == '#' % can have unlimited comments above the header
	line=fgetl(fid);
end
% parse the header line into # of traces, # of samples, dt, aligned phase, shift
dum=str2num(line); %#ok<ST2NM>
ntr=dum(1); nsamp=dum(2); dt=dum(3); align=dum(4); shift=dum(5);

% prep traces output for speed
traces = zeros(nsamp,3,ntr);
% Read each trace
for itr=1:ntr
    tr = textscan(fid,'%f %f %f',nsamp+2,'headerlines',3); % assume 3 lines of header separate each trace 
    traces(:,:,itr)= [tr{1},tr{2},tr{3}];
end
fclose(fid);

tt=(0:dt:dt*(nsamp-1)) - shift;

% 
% subplot(3,1,1)
% contplot(squeeze(traces(2,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
% title('SV')
% subplot(3,1,2)
% contplot(squeeze(traces(3,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
% title('SH');
% subplot(3,1,3)
% contplot(squeeze(traces(1,:,:))',[1:ntr]',dt,1,nsamp,1,ntr,1);
% title('P');
% xlabel('Trace number');
% ylabel('Time (s)');
