function [Ez, H_rho, H_phi] = getDielectricCylinderFieldUnderLineSource(radius, ...
                                                      cylinder, ...
                                                      background, ...
                                                      source_location, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    % [Ez, H_rho, H_phi] = getDielectricCylinderFieldUnderLineSource(radius, 
    %                                                cylinder,
    %                                                background,
    %                                                source_location,
    %                                                sensor_location,
    %                                                frequency)
    %
    % Calculate the field scattered by a dielectric cylinder centered
    % at the origin, due to an infinitely-long current source with a
    % unity current of 1A. The following convention of the coordinate
    % is assumed.
    %
    %                           y
    %                         ^
    %                         |
    %                         |
    %                         ---------> x
    %                        O
    %
    % 
    % Input:
    %
    % radius            scalar to denote the radius of
    %                   the cylinder (m)
    % cylinder          object of DielectricMaterial
    % background        object of DielectricMaterial
    % source_location   2x1 vector in the form of [x; y] (m)
    % sensor_location   2x1 vector in the form of [x; y] (m)
    % frequency         Nx1 vector (Hz)
    %
    % Output:
    %
    % Ez                Nx1 vector (V/m)
    % H_rho             Nx1 vector (V/m)
    % H_phi             Nx1 vector (V/m)
    %
    
    % $Author:: kzhu                                            $
    % $Rev:: 1614                                               $
    % $Date:: 2011-03-27 11:09:08 -0400 (Sun, 27 Mar 2011)      $

    p = inputParser;
    p.addRequired('radius', @(x)x>0);
    p.addRequired('cylinder',@isobject);
    p.addRequired('background', @isobject);
    p.addRequired('source_location', @isvector);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('verbose', 0, @(x)x==0||x==1);
    p.parse(radius, cylinder, background, source_location, sensor_location, frequency, varargin{:});
    if (p.Results.verbose)
        disp(p.Results);
    end
    
    omega          = 2*pi*frequency;
    EPS_O          = 8.854187817620389e-12;
    epsilon        = EPS_O*getPermittivity(background, frequency);
    [phi_s, rho_s] = cart2pol(source_location(1), source_location(2));
    [phi, rho]     = cart2pol(sensor_location(1), sensor_location(2));

    Ez    = zeros(length(frequency), 1);
    H_rho = zeros(length(frequency), 1);
    H_phi = zeros(length(frequency), 1);
    for iFreq = 1:length(frequency)
        k_1 = getElectromagneticWaveNumber(background, frequency(iFreq));
        k_2 = getElectromagneticWaveNumber(cylinder, frequency(iFreq));
        x_1 = k_1*radius;
        x_2 = k_2*radius;
        
        N_max = getN_max(radius, {cylinder}, ...
                                       background, ...
                                       frequency(iFreq));
        nu    = (-N_max:+N_max)';
        
        factor= -k_1.^2./(4.*omega.*epsilon);
        if (rho < rho_s) % inside source
            if (rho < radius) % inside cylinder
                dn_num =  k_1.*besselj(nu,x_1).*besselh_derivative(nu,2,x_1) ...
                         -k_1.*besselh(nu,2,x_1).*besselj_derivative(nu,x_1);
                dn_den =  k_1.*besselj(nu,x_2).*besselh_derivative(nu,2,x_1) ...
                         -k_2.*besselh(nu,2,x_1).*besselj_derivative(nu,x_2);
                dn     =  dn_num./dn_den;
                
                Ez(iFreq)    = factor.*sum(dn.*besselh(nu,2,k_1*rho_s).*besselj(nu,k_2*rho).*exp(j*nu*(phi-phi_s)));
                H_rho(iFreq) = 1/(4*j)/rho.*sum(dn.*besselh(nu,2,k_1*rho_s).*besselj(nu,k_2*rho).*(j*nu).*exp(j*nu*(phi-phi_s)));
                H_phi(iFreq) =-k_2/(4*j).*sum(dn.*besselh(nu,2,k_1*rho_s).*besselj(nu,k_2*rho).*exp(j*nu*(phi-phi_s)));
            else % outside cylinder
                cn_num =  k_2.*besselj(nu,x_1).*besselj_derivative(nu,x_2) ...
                         -k_1.*besselj(nu,x_2).*besselj_derivative(nu,x_1);
                cn_den =  k_1.*besselj(nu,x_2).*besselh_derivative(nu,2,x_1) ...
                         -k_2.*besselh(nu,2,x_1).*besselj_derivative(nu,x_2);
                cn     =  cn_num./cn_den;
                
                Ez(iFreq)    = factor.*sum(cn.*besselh(nu,2,k_1*rho_s).*besselh(nu,2,k_1*rho).*exp(j*nu*(phi-phi_s)));
                H_rho(iFreq) = 1/(4*j)/rho.*sum(cn.*besselh(nu,2,k_1*rho_s).*besselh(nu,2,k_1*rho).*(j*nu).*exp(j*nu*(phi-phi_s)));
                H_phi(iFreq) =-k_1/(4*j).*sum(cn.*besselh(nu,2,k_1*rho_s).*besselh_derivative(nu,2,k_1*rho).*exp(j*nu*(phi-phi_s)));
            end
        else % outside source radius
            cn_num =  k_2.*besselj(nu,x_1).*besselj_derivative(nu,x_2) ...
                     -k_1.*besselj(nu,x_2).*besselj_derivative(nu,x_1);
            cn_den =  k_1.*besselj(nu,x_2).*besselh_derivative(nu,2,x_1) ...
                     -k_2.*besselh(nu,2,x_1).*besselj_derivative(nu,x_2);
            cn     =  cn_num./cn_den;
            
            Ez(iFreq)    = factor.*sum(cn.*besselh(nu,2,k_1*rho).*besselh(nu,2,k_1*rho_s).*exp(j*nu*(phi-phi_s)));
            H_rho(iFreq) = 1/(4*j)/rho.*sum(cn.*besselh(nu,2,k_1*rho).*besselh(nu,2,k_1*rho_s).*(j*nu).*exp(j*nu*(phi-phi_s)));
            H_phi(iFreq) =-k_1/(4*j).*sum(cn.*besselh(nu,2,k_1*rho).*besselh_derivative(nu,2,k_1*rho_s).*exp(j*nu*(phi-phi_s)));
        end
    end
    


