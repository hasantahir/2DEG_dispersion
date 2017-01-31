function [epsilon_r, sigma_e] = getPermittivity(this, frequency)
% [epsilon_r, sigma_e] = getPermittivity(this)
%
% [epsilon_r, sigma_e] = getPermittivity(this, frequency)
%
% Input:
% 
% frequency     Nx1 vector (Hz)
%
    
% $Author:: kzhu                                      $
% $Rev:: 1487                                         $
% $Date: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $
    
    epsilon_r = this.epsilon_r;
    sigma_e   = this.sigma_e;
    
    
    