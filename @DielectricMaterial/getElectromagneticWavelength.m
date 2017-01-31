function wavelength = getElectromagneticWavelength(this, frequency)
% wavelength = getElectromagneticWaveLength(this, frequency)
%
% Input
%
% frequency     Nx1 vector (Hz)
%


% $Author:: kzhu                                          $
% $Rev:: 1487                                             $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011)    $

    wave_number = getElectromagneticWaveNumber(this, frequency);
    wavelength = 2*pi./real(wave_number);
    
    