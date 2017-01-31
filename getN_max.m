function N_max = getN_max(radius, sphere, background, frequency)
% N_max = getN_max(radius, sphere, background, frequency) implements
% [Wiscombe1980] to determine the number of bessel's function to be
% evaluated.
%
%
% Input:
%
% radius       N x 1 vector in the descending order
% sphere       N x 1 DielectricMaterial/DebyeMaterial cell array in the
%              same order as radius
% background   DielectricMaterial
% frequency    M x 1 vector
%
% Output:
%
% N_max        M x 1 vector
%
    
% $Author:: kzhu                                              $
% $Rev:: 1442                                                 $
% $Date:: 2011-02-06 03:08:52 -0500 (Sun, 06 Feb 2011)        $
    
    
    radius    = reshape(radius, 1, length(radius));
    sphere    = reshape(sphere, 1, length(sphere));
    frequency = reshape(frequency, length(frequency), 1);
    
    % Follow the index convention in [Yang2003]. The sphere is
    % specified from the inner most to the outter most.
    radius    = fliplr(radius);
    sphere    = fliplr(sphere);
    
    k_m       = conj(getElectromagneticWaveNumber(background, frequency));
    x         = k_m*radius; % x is a matrix
    x         = abs(x);     % x may be complex the background is lossy

    N_m       = conj(getComplexRefractiveIndex(background, frequency));    
    m         = zeros(length(frequency), length(sphere));
    for iSphere = 1:length(sphere)
        % Need to take the conjugate of the refractive index since I use
        % exp(j\omega t) as the time-harmonic factor while [Yang2003]
        % uses exp(-j\omega t).
        m(:, iSphere)  = conj(getComplexRefractiveIndex(sphere{iSphere}, frequency))./N_m;
    end

    
    N_stop      = ones(length(frequency), 1);
    idx         = find(0.02 <= x(:,end) & x(:,end) < 8);
    N_stop(idx) = round(x(idx,end)+4   *x(idx,end).^(1/3)+1); 
    idx         = find(8.00 <= x(:,end) & x(:,end) < 4200);
    N_stop(idx) = round(x(idx,end)+4.05*x(idx,end).^(1/3)+2);
    idx         = find(4200 <= x(:,end) & x(:,end) < 20000);
    N_stop(idx) = round(x(idx,end)+4   *x(idx,end).^(1/3)+2);

    N_max     = max([N_stop ...
                     round(abs(m.*x)) ...
                     round(abs(m(:,2:end).*x(:,1:end-1)))], [], 2)+15;

