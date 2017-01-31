function [Ez, H_rho, H_phi] = getCylindricalWaveUsingCylindricalExpansion(background, ...
                                                      source_location, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    % [Ez, H_rho, H_phi] = getCylindricalWaveUsingCylindricalExpansion(background,
    %                                                  source_location,
    %                                                  sensor_location, 
    %                                                  frequency) 
    %
    % Calculate the cylindrical wave due to an infinitely-long current
    % source with a unity current of 1A. This is the numerical
    % implementation of (5-119) based on (5-103) in [Harrington2001]
    %
    % 
    % Input:
    % 
    % background         DielectricMaterial
    % source_location    2x1 vector in the form of [x; y] (m)
    % sensor_location    2x1 vector in the form of [x; y] (m)
    % frequency          Nx1 vector in (Hz)
    %
    %
    % Output:
    %
    % Ez                  Nx1 vector (V/m)
    % H_rho               Nx1 vector (A/m)
    % H_phi               Nx1 vector (A/m)
    
    % $Author:: kzhu                                       $
    % $Rev:: 1613                                          $
    % $Date:: 2011-03-27 00:14:11 -0400 (Sun, 27 Mar 2011) $

    p = inputParser;
    p.addRequired('background', @isobject);
    p.addRequired('source_location', @isvector);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('verbose', 0, @(x)x==0||x==1);
    p.parse(background, source_location, sensor_location, frequency, varargin{:});
    if (p.Results.verbose)
        disp(p.Results);
    end

    omega          = 2*pi*frequency;
    EPS_O          = 8.854187817620389e-12;
    [phi_s, rho_s] = cart2pol(source_location(1), source_location(2));
    [phi, rho]     = cart2pol(sensor_location(1), sensor_location(2));
    
    Ez        = zeros(length(frequency), 1);
    H_rho     = zeros(length(frequency), 1);
    H_phi     = zeros(length(frequency), 1);
    for iFreq = 1:length(frequency);
        k     = getElectromagneticWaveNumber(background, ...
                                             frequency(iFreq)); % complex wavenumber of the background
        epsilon= EPS_O*getPermittivity(background, ...
                                       frequency(iFreq)); % complex wavenumber of the background
        N_max = getN_max(rho, {background}, ...
                              background, ...
                              frequency(iFreq));
        
        nu    = -N_max:+N_max;
        factor= -k.^2./(4.*omega.*epsilon);
        
        x     = k*rho;
        x_s   = k*rho_s;

        if (rho < rho_s) % inside source
            Ez(iFreq)    = factor.*sum(besselh(nu,2,x_s).*besselj(nu,x).*exp(j*nu*(phi-phi_s)));
            H_rho(iFreq) = 1/(4*j)/rho.*sum(besselh(nu,2,x_s).*besselj(nu,x).*(j*nu).*exp(j*nu*(phi-phi_s)));
            H_phi(iFreq) =-k/(4*j).*sum(besselh(nu,2,x_s).*besselj_derivative(nu,x).*exp(j*nu*(phi-phi_s)));
        else % outside source circle
            Ez(iFreq)    = factor.*sum(besselh(nu,2,x).*besselj(nu,x_s).*exp(j*nu*(phi-phi_s)));
            H_rho(iFreq) = 1/(4*j)/rho.*sum(besselh(nu,2,x).*besselj(nu,x_s).*(j*nu).*exp(j*nu*(phi-phi_s)));
            H_phi(iFreq) =-k/(4*j).*sum(besselh_derivative(nu,2,x).*besselj(nu,x_s).*exp(j*nu*(phi-phi_s)));
        end
    end
    



    
    
