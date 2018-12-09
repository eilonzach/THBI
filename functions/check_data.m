function check_data(trudata,par)
% check_data(trudata)
%
% function to check the data that was just loaded in, ascertaining whether
% or not it is usable and throwing up an error if not

if isempty(trudata)
    fprintf('No data for station %s, nwk %s\n\n',par.data.stadeets.sta,par.data.stadeets.nwk);
	error('No data for station %s, nwk %s\n\n',par.data.stadeets.sta,par.data.stadeets.nwk);
end

allpdytp = parse_dtype_all(par);
if ~any(strcmp({allpdytp{:,1}},'SW'))
    fprintf('No surface wave data at all for station %s, nwk %s\n\n',par.data.stadeets.sta,par.data.stadeets.nwk);
	error('No surface wave data at all for station %s, nwk %s\n\n',par.data.stadeets.sta,par.data.stadeets.nwk);
end


end

