function [traces,tt] = readPROPMATtr(trfile)
% [traces,tt] = readPROPMATtr(trfile)

% read one header to get nsamps etc.
fid=fopen(trfile,'r');
fgetl(fid); % skip one line for component name
line=fgetl(fid); % read line with nsamp data etc.
dum=str2num(line); %#ok<ST2NM>
nsamp=dum(1); dt=dum(2);
fclose(fid);
% lines of data per component
nl = nsamp/8;
% prep output structure
traces = nan(nsamp,3);
for ic = 1:3
    nh = 2 + (ic-1)*(nl+2); % number of header lines (inc previous components and their headers)
    fid=fopen(trfile,'r');
    C = textscan(fid,' %7d %7d %7d %7d %7d %7d %7d %7d',nl,'Headerlines',nh);
    fclose(fid);
    % convert cell to [nsamp/8 * 8] matrix, transpose and cat along columns
    try
        traces(:,ic) = reshape(cell2mat(C)',nsamp,1);
    catch
        traces(:,ic) = nan;
        fprintf('WARNING BAD PROPMAT OUTPUT\n')
    end
end



tt=(0:dt:dt*(nsamp-1))';

% 
% 
% figure(2);
% plot(tt,traces)

%% old version
% trs = cell(1,3);
% 
% %% VERSION 1
% fid=fopen(trfile,'r');
% % Read header:
% for ic = 1:3
% line=fgetl(fid); % dummy line for component name
% line=fgetl(fid);
% % parse the header line into "nsamps,dt,0,0,0"
% dum=str2num(line); %#ok<ST2NM>
% nsamp=dum(1); dt=dum(2);
% tr = zeros(nsamp,1);
% for il = 1:nsamp/8
%     line=fgetl(fid);
%     line = regexprep(line,'*******','NaN');
%     dum=str2num(line); %#ok<ST2NM>
%     tr([1:8] + 8*(il-1))=dum(:);
% end
% trs{ic} = tr;
% 
% end
% fclose(fid);
% traces = [trs{1},trs{2},trs{3}];
