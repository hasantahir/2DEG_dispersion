function [E_z, H_rho, H_phi] = getConductingCylinderFieldUnderPlaneWave(radius, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin);
    % [E_z, H_rho, H_phi] = getConductingCylinderFieldUnderPlaneWave(radius,
    %                                               background,
    %                                               sensor_location,
    %                                               frequency)
    %
    % Calculate the field of scattered by a conducting cylinder due to
    % an incident TM_z plane wave. It is the numerical implementation
    % of (5-107) in [Harrington2001].
    %
    % Input:
    %
    % radius           scalar to denote the radius
    %                  of the cylinder (m)
    % background       object of DielectricMaterial
    % sensor_location  2x1 vector in the form of [x; y] (m)
    % frequency        Nx1 vector (Hz)
    %
    % Output:
    %
    % E_z              Nx1 vector (V/m)
    % H_rho            Nx1 vector (A/m)
    % H_phi            Nx1 vector (A/m)


    % $Author:: kzhu                                       $
    % $Rev:: 1655                                          $
    % $Date:: 2011-04-13 15:44:49 -0400 (Wed, 13 Apr 2011) $

    nu    = -50:50; % order of the Bessel function
    k     = getElectromagneticWaveNumber(background, frequency);
    omega = 2*pi*frequency;
    mu    = 4*pi*1e-7*getPermeability(background);
    rho   = norm(sensor_location);
    if (rho < radius)
        E_z   = 0;
        H_rho = 0;
        H_phi = 0;
    else
        phi = angle(sensor_location(1)+j*sensor_location(2));
        x   = k*radius;
        a_n = -besselj(nu,x)./besselh(nu,2,x);
        E_z = sum(j.^(-nu).*a_n.*besselh(nu,2,k*rho).*exp(j*nu*phi));
        
        H_rho = -1./(j*omega*mu)/rho*sum(j.^(-nu)*j.*nu.*a_n.*besselh(nu,2,k*rho).*exp(j*nu*phi));
        H_phi = -1./(j*omega*mu)*(-1).*k*sum(j.^(-nu).*a_n.*besselh_derivative(nu,2,k*rho).*exp(j*nu*phi));
    end
