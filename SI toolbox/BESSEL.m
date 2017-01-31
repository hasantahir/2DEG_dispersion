function [J0,J1] = BESSEL(z)
%%
% Return the bessel function of the first kind and its derivative
% INPUT:
%       z - input complex number as the argument of the bessel function
%       J0 - cylindrical bessel function of the first kind
%       J1 - derivative of the cylindrical bessel function of the first
%       kind (J_0'(x) = -J_1(x))
% OUTPUT:
%       z_root - \sqrt(z) according to the vertical branch cut (xk)
%
% Hasan Tahir Abbas
% 01-14-2016
J0 =  besselj(0,z);
J1 =  -besselj(1,z);
end