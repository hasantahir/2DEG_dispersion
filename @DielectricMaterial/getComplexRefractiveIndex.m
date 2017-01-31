function ref_idx = getComplexRefractiveIndex(this, frequency)
% ref_idx = getComplexRefractiveIndex(this, frequency) computes the
% refractive index.
%
% Input:
%
% frequency    Nx1 vector (Hz)
%
    
% $Author:: kzhu                                      $
% $Rev:: 1487                                         $
% $Date: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $

    eps_r = getComplexPermittivity(this, frequency);
    mu_r  = getComplexPermeability(this, frequency);
    ref_idx =  sqrt(eps_r.*mu_r);
    

