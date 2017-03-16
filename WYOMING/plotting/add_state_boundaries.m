function add_state_boundaries(ax,latlims,lonlims,col)
% add_state_boundaries(ax,latlims,lonlims)
%   function to plot light grey state boundaries on maps

if isempty(ax)
    ax = gca;
end
if nargin<4 || isempty(col)
    col = [0.6 0.6 0.6];
end

latlims = latlims(:); % make column
lonlims = lonlims(:); % make column

states = shaperead('usastatehi','UseGeoCoords', true, 'BoundingBox', [lonlims,latlims]);

geoshow(ax, states, 'Edgecolor', col,'linestyle','-','facecolor','none')



end

