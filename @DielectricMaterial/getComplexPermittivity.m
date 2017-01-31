function epsilon_r = getComplexPermittivity(this, frequency)
% epsilon_r = getComplexPermittivity(this, frequency) computes the
% relative complex permittivity.
%
% Input:
%
%  frequency   Nx1 vector (Hz)
%
%
    
% $Author:: kzhu                                       $
% $Rev:: 1487                                          $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $
    
    
    EPS_O          = 8.8541878176e-12;    
    idx            = find(frequency==0);
    epsilon_r      = this.epsilon_r + this.sigma_e./(j*2*pi*frequency*EPS_O);
    epsilon_r(idx) = this.epsilon_r;
    
    
