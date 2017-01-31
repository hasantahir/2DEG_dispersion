function mu_r = getComplexPermeability(this, frequency)
% mu_r = getComplexPermeability(this, frequency) computes the relative
% complex permeability.
%
% Input:
%
% frequency    Nx1 vector (Hz)
%
    
     
% $Author:: kzhu                                        $
% $Rev:: 1487                                           $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011)  $

    
    MU_O      = 4*pi*1e-7;
    idx       = find(frequency==0);
    mu_r      = this.mu_r + this.sigma_m./(j*2*pi*frequency*MU_O);
    mu_r(idx) = this.mu_r;
    
    
    
    