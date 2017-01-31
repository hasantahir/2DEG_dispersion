% File name: visser.m
% Contains data for five-layer waveguide with gain and losses
% Reference:
% T.D. Visser et al, JQE v.31, p.1803 (1995)
% Fig.6
% Global variables to be transferred to function f_TE.m
%
global n_c % ref. index cladding
global n_layer % ref. index of internal layers
global n_s % ref. index substrate
global d_c % thickness of cladding (microns)
global d_layer % thicknesses of internal layers (microns)
global d_s % thickness of substrate (microns)
global k_0 % wavenumber
global NumberMesh % number of mesh points in each layer
% (including substrate and cladding)
n_c = 1.0;
n_layer = [3.40-1i*0.002 3.60+1i*0.010 3.40-1i*0.002];
n_s = 1.0;
d_s = 0.4;
d_layer = [0.6 0.4 0.6];
d_c = 0.5;
NumberMesh = [10 10 10 10 10];
lambda = 1.3; % wavelength in microns
k_0 = 2*pi/lambda;