function x = mesh_x(d_s,d_layer,d_c,NumberMesh)
% Generates one-dimensional mesh along x-axis
% Variable description:
% Input
% d_layer - thicknesses of each layer
% NumberMesh - number of mesh points in each layer
% Output
% x - mesh point coordinates
%
d_total = [d_s,d_layer,d_c]; % thicknesses of all layers
NumberOfLayers = length(d_total); % determine number of layers
delta = d_total./NumberMesh; % separation of points for all layers
%
x(1) = 0.0; % coordinate of first mesh point
i_mesh = 1;
for k = 1:NumberOfLayers % loop over all layers
    for i = 1:NumberMesh(k) % loop within layer
        x(i_mesh+1) = x(i_mesh) + delta(k);
        i_mesh = i_mesh + 1;
    end
end
