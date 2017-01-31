function H = besselh_derivative(nu,K,x)
% besselh_derivative(nu,K,x) using the recursive relationship to calculate the derivative of
% bessel's funtion of third kind, i.e. Hankel's function
%
% nu order of the bessel's function
% K = 1 if it is Hankel's function of the first kind; K=2 if it is Hankel's function of the
% second kind'
% x can be either real or complex
%
%
    
%   $Rev:: 580                                           $
%   $Author:: kzhu                                       $
%   $Date:: 2010-02-02 22:39:39 -0500 (Tue, 02 Feb 2010) $


    if (K ~=1 && K~=2)
        error('Improper kind of Hankel function');
    end
    
    H = 0.5*(besselh(nu-1,K,x)-besselh(nu+1,K,x));
end 


