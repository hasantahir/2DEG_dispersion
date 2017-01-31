function  [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                      sphere, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    
    % [E_r, E_theta, E_phi, H_r, H_theta, H_phi] =
    % getDielectricSphereFieldUnderPlaneWave(radius, ...
    %                                        sphere, ...
    %                                        background, ...
    %                                        sensor_location, ...
    %                                        frequency)
    %
    % Calculate the field scattered by a dielectric sphere centered at
    % the origine due to an incident x-polarized plane wave. The
    % scattered field is in (11-239) in [Balanis1989].  see the notes
    % on 2008-05-24 for the coefficients, a_n, b_n, and c_n.
    %
    %
    % See Fig. 11-25 in [Balanis1989] for the exact geometry.
    %    
    %
    % Input:
    %
    % radius             scalar to denote the radius of the sphere (m)
    % sphere             object of DielectricMaterial
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
    % $Rev:: 1704                                                  $
    % $Date:: 2011-04-17 20:19:47 -0400 (Sun, 17 Apr 2011)         $


    p = inputParser;    
    p.addRequired('radius', @isnumeric);
    p.addRequired('sphere', @isobject);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParamValue('debug', 0, @(x)x==0||x==1);
    p.parse(radius, sphere, background, sensor_location, frequency, varargin{:});
    if (p.Results.debug)
        disp(p.Results);
    end
    
    EPS_O     = 8.8541878176*1e-12;    
    MU_O      = 4*pi*1e-7;

    nFreq   = length(frequency);

    omega   = 2*pi*reshape(frequency,[nFreq 1]);    
    eta     = getIntrinsicImpedance(background, frequency);
    k       = getElectromagneticWaveNumber(background, frequency);
    mu      = getComplexPermeability(background, frequency)*MU_O;
    eps     = getComplexPermittivity(background, frequency)*EPS_O;
    eta_d   = getIntrinsicImpedance(sphere, frequency);
    k_d     = getElectromagneticWaveNumber(sphere, frequency);
    mu_d    = getComplexPermeability(sphere, frequency)*MU_O;
    eps_d   = getComplexPermittivity(sphere, frequency)*EPS_O;

    N       = getN_max(p.Results.radius, ...
                       {p.Results.sphere}, ...
                       p.Results.background, ...
                       p.Results.frequency);
    N       = max(N);
    nu      = 1:N;

    [r, theta, phi] = cart2sph(sensor_location(1),sensor_location(2),sensor_location(3));

    % Compute coefficients 
    a_n = j.^(-nu).*(2*nu+1)./(nu.*(nu+1)); 
    a_n = ones(nFreq,1)*a_n;

    % temp2 denotes the expression
    % kzlegendre(nu,1,cos(theta))/sin(theta). Here I am using a
    % recursive relation to compute temp2, which avoids the numerical
    % difficulty when theta == 0 or PI.
    temp2    = zeros(1, length(nu)); 
    temp2(1) = -1; 
    temp2(2) = -3*cos(theta);  
    for n = 2:nu(end-1) 
        temp2(n+1) = (2*n+1)/n*cos(theta)*temp2(n) - (n+1)/n*temp2(n-1);
    end
    
    % temp1 denotes the expression
    % sin(theta)*kzlegendre_derivative(nu,1,cos(theta)).  Here I am
    % also using a recursive relation to compute temp1 from temp2,
    % which avoids numerical difficulty when theta == 0 or PI.
    temp1    = zeros(1, length(nu));
    temp1(1) = cos(theta); 
    for n = 2:nu(end)
        temp1(n) = (n+1)*temp2(n-1)-n*cos(theta)*temp2(n);
    end

    temp1     = ones(nFreq,1)*temp1;
    temp2     = ones(nFreq,1)*temp2;

    iNU = 10;    
    if (p.Results.debug);
        
        A   = [ric_besselh_derivative(iNU,2,k*radius) -sqrt(eps*mu)/sqrt(eps_d*mu_d)*ric_besselj_derivative(iNU,k_d*radius);
               ric_besselh(iNU,2,k*radius) -mu/mu_d*ric_besselj(iNU,k_d*radius)];
        rhs = -a_n(iNU)*[ric_besselj_derivative(iNU,k*radius); ric_besselj(iNU,k*radius)];
        x   = A\rhs;
        disp(['b_n ' num2str(x(1)) ' d_n ' num2str(x(2))]);

        A   = [ric_besselh(iNU,2,k*radius) -sqrt(eps*mu)/sqrt(eps_d*mu_d)*ric_besselj(iNU,k_d*radius);
               ric_besselh_derivative(iNU,2,k*radius) -mu/mu_d*ric_besselj_derivative(iNU,k_d*radius)];
        rhs = -a_n(iNU)*[ric_besselj(iNU,k*radius); ric_besselj_derivative(iNU,k*radius)];
        x   = A\rhs;
        disp(['c_n ' num2str(x(1)) ' e_n ' num2str(x(2))]);
        disp('------');
    end
    
    
    if (r < radius)
        num = j.*mu_d/sqrt(mu)*sqrt(eps_d);
        den =  - sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselh_derivative(nu,2,k*radius))...
               + sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselh(nu,2,k*radius)).*transpose(ric_besselj_derivative(nu,k_d*radius));
        d_n = num*ones(1,N)./den.*a_n;
        
        den = + sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselh(nu,2,k*radius)).*transpose(ric_besselj_derivative(nu,k_d*radius))...
              - sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselh_derivative(nu,2,k*radius));
        e_n = num*ones(1,N)./den.*a_n;
        

        if (p.Results.debug)
            disp(['d_n ' num2str(d_n(iNU)) ' e_n ' num2str(e_n(iNU))]);
            return
        end
        
        x   = k_d*r;
        
        % Implement (11-239a) in [Balanis1989]    
        alpha     = (transpose(ric_besselj_derivative(nu,x,2))+transpose(ric_besselj(nu,x)))...
            .*transpose(kzlegendre(nu,1,cos(theta))*ones(1,nFreq));
        E_r       = -j*cos(phi)*sum(d_n.*alpha, 2);
        H_r       = -j*sin(phi)*sum(e_n.*alpha, 2)./eta_d;
        
        % Implement (11-239b) in [Balanis1989]    
        alpha     = transpose(ric_besselj_derivative(nu,x)).*temp1;
        beta      = transpose(ric_besselj(nu,x)).*temp2;
        summation = j*d_n.*alpha - e_n.*beta;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = j*e_n.*alpha - d_n.*beta;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta_d;
        
        % Implement (11-239c) in [Balanis1989]
        alpha     = transpose(ric_besselj_derivative(nu,x)).*temp2;
        beta      = transpose(ric_besselj(nu,x)).*temp1;
        summation = j*d_n.*alpha - e_n.*beta;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = j*e_n.*alpha - d_n.*beta;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta_d;        
    else
        
        num =  + sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselj(nu,k*radius))  .*transpose(ric_besselj_derivative(nu,k_d*radius)) ...
               - sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselj_derivative(nu,k*radius));
        den =  + sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselh_derivative(nu,2,k*radius))...
               - sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselh(nu,2,k*radius)).*transpose(ric_besselj_derivative(nu,k_d*radius));
        b_n = num./den.*a_n;
        
        num = + sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselj_derivative(nu,k*radius))...
              - sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselj(nu,k*radius))  .*transpose(ric_besselj_derivative(nu,k_d*radius));
        den = + sqrt(mu.*eps_d)*ones(1,N).*transpose(ric_besselh(nu,2,k*radius)).*transpose(ric_besselj_derivative(nu,k_d*radius))...
              - sqrt(mu_d.*eps)*ones(1,N).*transpose(ric_besselj(nu,k_d*radius)).*transpose(ric_besselh_derivative(nu,2,k*radius));
        c_n = num./den.*a_n;
        
        if (p.Results.debug)
            disp(['b_n ' num2str(b_n(iNU)) ' c_n ' num2str(c_n(iNU))]);
            return
        end
        
        x   = k*r;
        
        % Implement (11-239a) in [Balanis1989]    
        alpha     = (transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
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
    end