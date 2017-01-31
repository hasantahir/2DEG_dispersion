function [mu_r, sigma_m] = getPermeability(this, frequency)
% [mu_r, sigma_m] = getPermeability(this, frequency)
%
% Input:
% 
% frequency    Nx1 vector (Hz)
    
    
% $Author:: kzhu                                       $
% $Rev:: 1487                                          $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $
    
    mu_r    = this.mu_r;
    sigma_m = this.sigma_m;
    