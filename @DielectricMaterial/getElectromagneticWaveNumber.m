function wave_number = getElectromagneticWaveNumber(this, frequency)
% wave_number = getElectromagneticWaveNumber(this, frequency)
%
% Input:
% 
% frequency     Nx1 vector (Hz)
%
%
    
% $Author:: kzhu                                          $
% $Rev:: 1487                                             $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011)    $
    

    C_O =      2.997924580003452e+08;
    permittivity = getComplexPermittivity(this, frequency);
    permeability = getComplexPermeability(this, frequency);
    wave_number  = 2*pi*frequency.*sqrt(permittivity.*permeability)/C_O;
    
    