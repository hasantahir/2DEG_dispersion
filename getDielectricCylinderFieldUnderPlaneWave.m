function [Ez, H_rho, H_phi] = getDielectricCylinderFieldUnderPlaneWave(radius, ...
                                                      cylinder, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    % [Ez, H_rho, H_phi] = getDielectricCylinderFieldUnderPlaneWave(radius, cylinder,
    %                                               background,
    %                                               sensor_location,
    %                                               frequency)
    %
    % Calculate the field of scattered by a dielectric cylinder due to
    % an incident TM_z plane wave. If the sensor_location is outside
    % the cylinder, this function returns the scattered field. If the
    % sensor_location is inside the cylinder, it returns the total
    % field inside the cylinder. The incident TM_z plane wave is
    % assumed propagating from left (neg X) to right (pos X).
    %
    % This problem is solved in Problem 11-26 in [Balanis1989]. 
    %
    % Input:
    %
    % radius            scalar to denote the radius of
    %                   the cylinder (m)
    % cylinder          object of DielectricMaterial
    % background        object of DielectricMaterial
    % sensor_location   2x1 vector in the form of [x; y] (m)
    % frequency         Nx1 vector (Hz)
    %
    % Output:
    %
    % Ez                Nx1 vector (V/m)
    % H_rho             Nx1 vector (A/m)
    % H_phi             Nx1 vector (A/m)
    % 
    
    % $Author:: kzhu                                       $
    % $Rev:: 1611                                          $
    % $Date:: 2011-03-26 11:10:33 -0400 (Sat, 26 Mar 2011) $


    p = inputParser;
    p.addRequired('radius', @(x)x>0);
    p.addRequired('cylinder',@isobject);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('verbose', 0, @(x)x==0||x==1);
    p.parse(radius, cylinder, background, sensor_location, frequency, varargin{:});
    if (p.Results.verbose)
        disp(p.Results);
    end
    
    [phi, rho]= cart2pol(p.Results.sensor_location(1), p.Results.sensor_location(2));


    Ez        = zeros(length(p.Results.frequency), 1);
    H_rho     = zeros(length(p.Results.frequency), 1);
    H_phi     = zeros(length(p.Results.frequency), 1);
    for iFreq = 1:length(p.Results.frequency)
        omega = 2*pi*frequency;
        k     = getElectromagneticWaveNumber(p.Results.background, ...
                                             p.Results.frequency(iFreq));
        eta   = getIntrinsicImpedance(p.Results.background, ...
                                      p.Results.frequency(iFreq));
        k_d   = getElectromagneticWaveNumber(p.Results.cylinder, ...
                                             p.Results.frequency(iFreq));
        eta_d = getIntrinsicImpedance(p.Results.cylinder, ...
                                      p.Results.frequency(iFreq));
        
        N_max = getN_max(p.Results.radius, ...
                         {p.Results.cylinder}, ...
                         p.Results.background, ...
                         p.Results.frequency(iFreq));
        n     = (-N_max:+N_max)';
        
        x = k_d*radius;
        y = k*radius;
        
        temp  = (min(n)-1:max(n+1))';
        a = besselj(temp,x); d = 0.5 * (a(1:end-2)-a(3:end));
        b = besselj(temp,y); e = 0.5 * (b(1:end-2)-b(3:end));
        c = besselh(temp,2,y);f= 0.5 * (c(1:end-2)-c(3:end));
        a = a(2:end-1);
        b = b(2:end-1);
        c = c(2:end-1);
        
        if (rho > radius) % scatter
            num = 1/eta_d.*d.*b-1/eta.*e.*a;
            den = 1/eta.*a.*f -1/eta_d.*d.*c;
            a_n = num./den.*j.^(-n);
            
            Ez(iFreq)    =              sum(a_n.*besselh(n,2,k*rho).*exp(j*n*phi));
            H_rho(iFreq) = -1./(j*eta).*sum(a_n.*besselh(n,2,k*rho)./(k*rho).*j.*n.*exp(j*n*phi));
            H_phi(iFreq) = +1./(j*eta).*sum(a_n.*besselh_derivative(n,2,k*rho).*exp(j*n*phi));
        else % inside
            num = 1/eta.*f.*b-1/eta.*e.*c;
            den = 1/eta.*f.*a-1/eta_d.*d.*c;
            c_n = num./den.*j.^(-n);

            Ez(iFreq)    =                sum(c_n.*besselj(n,k_d*rho).*exp(j*n*phi));
            H_rho(iFreq) = -1./(j*eta_d).*sum(c_n.*besselj(n,k_d*rho)./(k_d*rho).*j.*n.*exp(j*n*phi));
            H_phi(iFreq) = +1./(j*eta_d).*sum(c_n.*besselj_derivative(n,k_d*rho).*exp(j*n*phi));
        end
    end


