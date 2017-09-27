function [ model ] = expand_dathparms_to_data( model,trudata,par )
%[ model ] = expand_dathparms_to_data( model,trudata )
%   Function to expand the fields of the hyperparameters in the model to 
%   allow for multiple traces in each datatype.

fni = par.inv.datatypes;
fnd = fieldnames(trudata);
fnm = fieldnames(model.datahparm);

% cycle through data, removing datatypes if not in par.
for jfnd = 1:length(fnd)
    dtype = regexprep(fnd{jfnd},'sig_','');
    if ~any(strcmp(dtype,fni))
        trudata = rmfield(trudata,dtype);
    end
end
fnd = fieldnames(trudata);

% cycle through data, expanding datahparm fields where needed
for jfnd = 1:length(fnd)
    model.datahparm.(['sig_',fnd{jfnd}]) = ...
        ones(length(trudata.(fnd{jfnd})),1).*model.datahparm.(['sig_',fnd{jfnd}]);
end

% cycle through model, removing datatype sigmas if they aren't in the real
% data
for jfnm = 1:length(fnm)
    dtype = regexprep(fnm{jfnm},'sig_','');
    if ~any(strcmp(dtype,fnd))
        model.datahparm = rmfield(model.datahparm,fnm{jfnm});
    end
end



end

