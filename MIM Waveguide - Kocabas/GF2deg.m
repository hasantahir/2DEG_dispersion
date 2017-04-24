function D = GF2deg(kp)
% ***********************************************************************
%
%      Computes the Transmission Line Greens Function of a three
%      layer structure with source located at the center of the middle
%      layer
%
% ***********************************************************************


%% Calls polee_search_MIM
% Data from datafile of test_case3 from the reference paper
f = 25e12;
c = 3e8;
lambda = c/f;
omega = 2*pi*c/lambda;

% Example Validations

% Material Properties
ep1 = 8 + 1i*.36;
% ep2 = -0.0884 + 1i*.0005;
ep3 = -9.884 + 1i*1.9661;
cond =  1i*1e-5;



% EM constants
mu0 = 4*pi*1e-7;
ep0 = 8.854e-12;

% Propagations Constants
k1 = omega*sqrt(mu0*ep0*ep1);
k3 = omega*sqrt(mu0*ep0*ep3);


d = 0;






kz1 =  sqrtbr(k1 ^2  - kp .^2, pi/2);

kz3 =  -sqrtbr(k3 ^2  - kp .^2, pi/2);

% TE/TM switch
nu = 1;



% if imag(kz1) >= 0
%     kz1 = conj(kz1);
% end
% 
% 
% if imag(kz3) <= 0
%     kz3 = conj(kz3);
% end



    % TM case
    Z1 =  kz1./(omega*ep1*ep0);
    Z2 =  1./cond;
    Z3 =  kz3./(omega*ep3*ep0);





%% Reflection Coefficients
Gamma_left = (Z3 - Z2) ./ (Z3 + Z2); % Left-looking
Gamma_right =  (Z1 - Z2) ./ (Z1 + Z2); % Right-looking

% Denominator, Dispersion relation
% D =  1 - Gamma_left.*Gamma_right;
D = (Z1 + Z3)./(Z1.*Z3) - 1./Z2;

end