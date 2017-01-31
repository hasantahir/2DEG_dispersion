function [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getConductingSphereFieldUnderPlaneWave(radius, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    
    % [E_r, E_theta, E_phi, H_r, H_theta, H_phi] =
    % getConductingSphereFieldUnderPlaneWaveExcitation(radius, ...
    %                                                  background, ...
    %                                                  sensor_location, ...
    %                                                  background, ...
    %                                                  frequency)
    %
    % Calculate the field scattered by a conducting sphere centered at
    % the origine due to an incident x-polarized plane wave.  The
    % coefficients, a_n, b_n, and c_n are in (11-231a), (11-238a), and
    % (11-238b) in [Balanis1989].  And the scattered field is in
    % (11-239).
    %
    % See Fig. 11-25 in [Balanis1989] for the exact geometry.
    %    
    %
    % Input:
    %
    % radius             scalar to denote the radius of the sphere (m)
    % background         object of DielectricMaterial
    % sensor_location    3x1 vector in the form of [x; y; z] (m)
    % frequency          Nx1 vector in (Hz)
    % 
    % Output:
    %
    % E_r               Nx1 vector (V/m)
    % E_phi             Nx1 vector (V/m)
    % E_theta           Nx1 vector (V/m)
    % H_r               Nx1 vector (A/m)
    % H_phi             Nx1 vector (A/m)
    % H_theta           Nx1 vector (A/m)

    % $Author:: kzhu                                               $
    % $Rev:: 1696                                                  $
    % $Date:: 2011-04-17 19:02:10 -0400 (Sun, 17 Apr 2011)         $


    p = inputParser;    
    p.addRequired('radius', @isnumeric);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('debug', 0, @(x)x==0||x==1);
    p.parse(radius, background, sensor_location, frequency, varargin{:});
    if (p.Results.debug)
        disp(p.Results);
    end

    nFreq   = length(frequency);
    omega   = 2*pi*reshape(frequency,[nFreq 1]);
    eta     = getIntrinsicImpedance(background, frequency);
    k       = getElectromagneticWaveNumber(background, frequency);

    N       = 50; % Order of the Bessel(Hankel). Important. If this number is
                  % not large enough, the summing serioues does not converge.
    nu      = 1:N;

    [r, theta, phi] = cart2sph(sensor_location(1),sensor_location(2),sensor_location(3));

    if (r < radius)
        E_r = 0;
        E_theta = 0;
        E_phi = 0;
        H_r = 0;
        H_theta = 0;
        H_phi = 0;
        return;
    end

    % Compute coefficients as in (11-231a), (11-238a), (11-239b)
    a_n = j.^(-nu).*(2*nu+1)./(nu.*(nu+1)); 
    a_n = ones(nFreq,1)*a_n;
    b_n = -a_n.*transpose(ric_besselj_derivative(nu,k*radius))./transpose(ric_besselh_derivative(nu,2,k*radius));
    c_n = -a_n.*transpose(ric_besselj(nu,k*radius))./transpose(ric_besselh(nu,2, k*radius));
    
    % temp2 denote the expression kzlegendre(nu,1,cos(theta))/sin(theta). Here I
    % am using an recursive relation to compute temp2, which avoids the
    % numerical difficulty when theta == 0 or PI.
    temp2    = zeros(1, length(nu)); 
    temp2(1) = -1; 
    temp2(2) = -3*cos(theta);  
    for n = 2:nu(end-1) 
        temp2(n+1) = (2*n+1)/n*cos(theta)*temp2(n) - (n+1)/n*temp2(n-1);
    end
    
    % temp1 denote the expression
    % sin(theta)*kzlegendre_derivative(nu,1,cos(theta)).  Here I am
    % also using an recursive relation to compute temp1 from temp2,
    % which avoids numerical difficulty when theta == 0 or PI.
    temp1    = zeros(1, length(nu));
    temp1(1) = cos(theta); 
    for n = 2:nu(end)
        temp1(n) = (n+1)*temp2(n-1)-n*cos(theta)*temp2(n);
    end

    temp1     = ones(nFreq,1)*temp1;
    temp2     = ones(nFreq,1)*temp2;
    
    x   = k*r;
    
    % Implement (11-239a) in [Balanis1989]    
    alpha     = ( transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
        .*transpose(kzlegendre(nu,1,cos(theta))*ones(1,nFreq));
    E_r       = -j*cos(phi)*sum(b_n.*alpha, 2);
    H_r       = -j*sin(phi)*sum(c_n.*alpha, 2)./eta;
    
    % Implement (11-239b) in [Balanis1989]    
    alpha     = transpose(ric_besselh_derivative(nu,2,x)).*temp1;
    beta      = transpose(ric_besselh(nu,2,x)).*temp2;
    summation = j*b_n.*alpha - c_n.*beta;
    E_theta   = cos(phi)./x.*sum(summation,2);
    summation = j*c_n.*alpha - b_n.*beta;
    H_theta   = sin(phi)./x.*sum(summation,2)./eta;
    
    % Implement (11-239c) in [Balanis1989]
    alpha     = transpose(ric_besselh_derivative(nu,2,x)).*temp2;
    beta      = transpose(ric_besselh(nu,2,x)).*temp1;
    summation = j*b_n.*alpha - c_n.*beta;
    E_phi     = sin(phi)./x.*sum(summation,2);
    summation = j*c_n.*alpha - b_n.*beta;
    H_phi     =-cos(phi)./x.*sum(summation,2)./eta;


