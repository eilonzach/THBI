function trudata = duplicate_data_distribute(trudata,par)
% trudata = duplicate_data_distribute(trudata,par)
% distribute data for different processing (e.g. _lo, _cms)
    for idt = 1:length(par.inv.datatypes)
        dtype = par.inv.datatypes{idt}; 
        pdt = parse_dtype( dtype ); 
%         if ~isfield(trudata,par.inv.datatypes{idt}) && strcmp(pdt{1},'BW')
%             trudata.(dtype) = trudata.([pdt{1},'_',pdt{2}]); % insert standard BW if needed
%         end
        if  strcmp(pdt{4},'lo') && (strcmp(pdt{1},'BW') || strcmp(pdt{1},'RF'))
            if isfield(trudata,['BW_',pdt{2}])
                trudata.(regexprep(dtype,'RF','BW')) = trudata.(['BW_',pdt{2}]); % insert standard BW if needed
            end
        end
    end
end

