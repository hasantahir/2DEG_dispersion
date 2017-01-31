function [Ez, H_rho, H_phi] = getConductingCylinderFieldUnderLineSource(radius, ...
                                                      background, ...
                                                      source_location, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin) 
    % [Ez, H_rho, H_phi] = getConductingCylinderFieldUnderLineSource(radius,
    %                                                background,
    %                                                source_location,
    %                                                sensor_location,
    %                                                frequency)
    %
    % Calculate the field scattered by a conducting cylinder centered
    % at the origin, due to an infinitely-long current source with a
    % unity current of 1A. It is the numerical implementation of
    % (11-164) in [Balanis1989].
    %
    %
    % Input:
    % 
    % radius            scalar to denote the radius of
    %                   the cylinder (m)
    % background        object of DielectricMaterial
    % source_location   2x1 vector in the form of [x; y] (m)
    % sensor_location   2x1 vector in the form of [x; y] (m)
    % frequency         Nx1 vector (Hz)
    %
    % Output:
    %
    % Ez                Nx1 vector
    % H_rho             Nx1 vector
    % H_phi             Nx1 vector
    
    % $Author:: kzhu                                       $
    % $Rev:: 1614                                          $
    % $Date:: 2011-03-27 11:09:08 -0400 (Sun, 27 Mar 2011) $

    omega          = 2*pi*frequency;
    EPS_O          = 8.854187817620389e-12;
    epsilon        = EPS_O*getPermittivity(background);
    k              = getElectromagneticWaveNumber(background, frequency);
    [phi_s, rho_s] = cart2pol(source_location(1), source_location(2));
    [phi, rho]     = cart2pol(sensor_location(1), sensor_location(2));

    nu    = -50:50; 
    c_n   = -besselj(nu, k*radius)./besselh(nu,2,k*radius);

    factor= -k.^2./(4.*omega.*epsilon);
    if (rho < rho_s) % inside source
        if (rho < radius) % inside cylinder
            Ez    = 0;
            H_rho = 0;
            H_phi = 0;
        else % outside cylinder
            x     = k*rho;
            x_s   = k*rho_s;
            Ez    = factor.*sum(besselh(nu,2,x_s).*(besselj(nu,x)+c_n.*besselh(nu,2,x)).*exp(j*nu*(phi-phi_s)));
            H_rho = 1/(4*j)/rho.*sum(besselh(nu,2,x_s).*(besselj(nu,x)+c_n.*besselh(nu,2,x)).*(j*nu).*exp(j*nu*(phi-phi_s)));
            H_phi =-k/(4*j).*sum(besselh(nu,2,x_s).*(besselj_derivative(nu,x)+c_n.*besselh_derivative(nu,2,x)).*exp(j*nu*(phi-phi_s)));
        end
    else % outside source
        x     = k*rho;
        x_s   = k*rho_s;
        Ez    = factor.*sum(besselh(nu,2,x).*(besselj(nu,x_s)+c_n.*besselh(nu,2,x_s)).*exp(j*nu*(phi-phi_s)));
        H_rho = 1/(4*j)/rho.*sum(besselh(nu,2,x).*(besselj(nu,x_s)+c_n.*besselh(nu,2,x_s)).*(j*nu).*exp(j*nu*(phi-phi_s)));
        H_phi =-k/(4*j).*sum(besselh_derivative(nu,2,x).*(besselj(nu,x_s)+c_n.*besselh(nu,2,x_s)).*exp(j*nu*(phi-phi_s)));
    end
    
    
    
