function z_root = CROOT(z,xk,lambda)
%%
% INPUT:
%       z - input complex number to find correct square root
%       xk- the point in the complex lambda plane where the vertical branch
%           cut beings
%       lambda 
% OUTPUT:
%       z_root - \sqrt(z) according to the vertical branch cut (xk)
%
% Hasan Tahir Abbas
% 01-14-2016
z_root = sqrt(z);
if real(lambda) < real(xk) && imag(lambda) > imag(xk)
    z_root = - z_root;
end
end