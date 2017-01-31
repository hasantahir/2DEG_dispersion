function result = f_TE(z)
% Creates function used to determine propagation constant
% Variable description:
% result - expression used in search for propagation constant
% z - actual value of propagation constant
%
% Global variables:
% Global variables are used to transfer values from data functions
global n_s % ref. index substrate
global n_c % ref. index cladding
global n_layer % ref. index of internal layers
global d_layer % thicknesses of internal layers (microns)
global k_0 % wavenumber
%
zz=z*k_0;
NumLayers = length(d_layer);
%
% Creation for substrate and cladding
gamma_sub=sqrt(zz^2-(k_0*n_s)^2);
gamma_clad=sqrt(zz^2-(k_0*n_c)^2);
%
% Creation of kappa for internal layers
kappa=sqrt(k_0^2*n_layer.^2-zz.^2);
temp = kappa.*d_layer;
%
% Construction of transfer matrix for first layer
cc = cos(temp);
ss = sin(temp);
m(1,1) = cc(1);
m(1,2) = -1j*ss(1)/kappa(1);
m(2,1) = -1j*kappa(1)*ss(1);
m(2,2) = cc(1);
%
% Construction of transfer matrices for remaining layers
% and multiplication of matrices
for i=2:NumLayers
    mt(1,1) = cc(i);
    mt(1,2) = -1j*ss(i)/kappa(i);
    mt(2,1) = -1j*ss(i)*kappa(i);
    mt(2,2) = cc(i);
    m = mt*m;
end
%
result = 1j*(gamma_clad*m(1,1)+gamma_sub*m(2,2))...
    + m(2,1) - gamma_sub*gamma_clad*m(1,2);