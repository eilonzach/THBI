function allpdytp = parse_dtype_all(par)
% allpdytp = parse_dtype_all(par)
%    quickly parse all data types and put into structure that contains all
%    of them...

allpdytp = cell(length(par.inv.datatypes),4);
for id = 1:length(par.inv.datatypes)
    allpdytp(id,:)=parse_dtype(par.inv.datatypes{id});
end

end

