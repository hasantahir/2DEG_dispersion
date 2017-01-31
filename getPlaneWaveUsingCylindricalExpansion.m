function [Ex, Hy] = getPlaneWaveUsingCylindricalExpansion(background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    % [Ex, Hy] = getPlaneWaveUsingCylindricalExpasion(background,
    %                                                 sensor_location, 
    %                                                 frequency) 
    %
    % Calculate a plane wave using (5-101) in [Harrington2001]
    % assuming it propagating in the +z direction polarized in the +x
    % direction.
    % 
    % Input:
    % 
    % background          DielectricMaterial
    % sensor_location     [x;y;z] (m)
    % frequency           Nx1 vector (Hz)
    %
    % Optional:
    %
    % 'time'              scalar (s). Computes the field at a
    %                     particular time instance. Default to [].
    %
    % Output:
    %
    % Ex                  Nx1 vector (V/m)
    % Hy                  Nx1 vector (A/m)
    
    % $Author:: kzhu                                       $
    % $Rev:: 1696                                          $
    % $Date:: 2011-04-17 19:02:10 -0400 (Sun, 17 Apr 2011) $

    p = inputParser;
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('time', [], @isnumeric);
    p.parse(background, sensor_location, frequency, varargin{:});


    frequency = reshape(frequency, length(frequency), 1);
    [phi, rho]= cart2pol(sensor_location(3), sensor_location(2));

    k_m   = getElectromagneticWaveNumber(background, frequency);
    N_max = getN_max(rho, {background}, background, frequency);
    N_max = max(N_max);
    n     = (-N_max:+N_max);
    temp  = ones(length(frequency),1)*(j.^(-n))   ...
            .*besselj(n, k_m.*rho) ...
            .*(ones(length(frequency),1)*exp(j*n*phi));
    Ex    = sum(temp,2);
    Hy    = Ex./getIntrinsicImpedance(background, frequency);

    if (~isempty(p.Results.time))
        omega = 2*pi*frequency;
        Ex    = real(Ex*exp(j*omega*p.Results.time));
        Hy    = real(Hy*exp(j*omega*p.Results.time));
    end
    
    
    