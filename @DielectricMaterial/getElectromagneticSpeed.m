function prop_speed = getElectromagneticSpeed(this, frequency)
% prop_speed = getElectromagneticSpeed(this, frequency)
%
% Input:
%
% frequency    Nx1 vector (Hz)

% $Author:: kzhu                                       $
% $Rev:: 1569                                          $
% $Date:: 2011-03-16 17:50:52 -0400 (Wed, 16 Mar 2011) $
    
    C_O          = 2.997924580003452e+08;
    switch  (nargin)
      case 1
        prop_speed = C_O/sqrt(this.epsilon_r*this.mu_r);
      otherwise
        permittivity = getComplexPermittivity(this, frequency);
        permeability = getComplexPermeability(this, frequency);
        prop_speed   = C_O./sqrt(permittivity.*permeability);
        prop_speed   = real(prop_speed);
    end    
    
    
    