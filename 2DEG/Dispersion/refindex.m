function n_mesh = refindex(x,interface,index_layer)
% Assigns the values of refractive indices to mesh points
% in all layers
% Input
% x() -- mesh points coordinates
% interface(n) -- number of mesh points in layers
% index_layer() -- refrective index in layers
% Output
% index_mesh -- refractive index for each mesh point
%
% Within a given layer, refractive index is assigned the same value.
% Loop scans over all mesh points.
% For all mesh points selected for a given layer, the same
% value of refractive index is assigned.
%
N_mesh = length(x);
NumberOfLayers = length(index_layer);
%
i_mesh = 1;
for k = 1:NumberOfLayers % loop over all layers
    for i = 1:interface(k) % loop within layer
        n_mesh(i_mesh+1) = index_layer(k);
        i_mesh = i_mesh + 1;
    end
end