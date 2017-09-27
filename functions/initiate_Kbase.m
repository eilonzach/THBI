function Kbase = initiate_Kbase
% Kbase = initiate_Kbase
% 
%  Function to initiate the base kernel structure, with all the
%  possiblities for data types (Ray/Love, ph/gr)

kstruc = struct('periods',[],'phV',[],'grV',[],'Kph',{},'Kgr',{});

Kbase = struct('modelk',[],'Ray',kstruc,'Lov',kstruc);

end

