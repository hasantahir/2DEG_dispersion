function J = besselj_derivative(nu,x)
% J = besselj_derivative(nu,K,x) using the recursive relationship to
% calculate the derivative of bessel's funtion of first kind
%
% nu    order of the bessel's function
% x     can be either real or complex
%
%
    
%   $Rev:: 580                                           $
%   $Author:: kzhu                                       $
%   $Date:: 2010-02-02 22:39:39 -0500 (Tue, 02 Feb 2010) $
    
    J = 0.5*(besselj(nu-1,x)-besselj(nu+1,x));
end %besselj_derivative
