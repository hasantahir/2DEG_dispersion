%% solution for dispersion relation for 4 layer 2deg structure

clear all

f = linspace(.1e12, 1e12, 1e1); % Frequency in Hz
% f = .5e12;
% f = 8.3e12; % Frequency
% f = 25e12;
%% Constants
e0 = 8.85e-12; % Free-space permittivity
c = 3e8; % speed of light
e = 1.60218e-19; % electron charge
me = 9.109e-31; % Electron mass

% free space wavenumber
w = 2*pi*f;
k0 = 2*pi*f/c;
f0 = 1e12;
k_0 = 2*pi*f0/c;
% layer material properties
e1 = 1; % Air
e2 = 11.88; % AlGaAs
e3 = 13.01; % GaAs

% dispersive from palik @ .5 THz
% e2 = 11.88; % AlGaAs
% e3 = 13.01; % GaAs
% dispersive from palik @ 8.3 THz
% e2 = 7.44 + 1i*.703; % AlGaAs
% e3 = -15.62 - 1i*4.561;% + 1i*.703; % GaAs

% dispersive from palik @ 25 THz
% e2 = 9.179 - 1i*.0496; % AlGaAs
% e3 = 10.79; % GaAs
% Layer geometry
d1 = 20e-9;
d2 = 50e-9;

% 2DEG Material properties
Ns = 4e12 *1e4;
ms = .067*me; % Effective mass
mu = 1e6*1e-4;
tau = 1.14e-10;

% Transverse propagation constants
% kz1 = @(kx) sqrt(k0^2 - kx.^2);
% kz2 = @(kx) sqrt(k0^2*e2 - kx.^2);
% kz3 = @(kx) sqrt(k0^2*e3 - kx.^2);
root = [];
sigma = Ns*e^2*tau/ms./(1 + 1i*w*tau);

for i = 1 : length(w)
    
    kz1 = @(kx) sqrt(kx.^2 + k0(i)^2);
    kz2 = @(kx) sqrt(kx.^2 + k0(i)^2*e2);
    kz3 = @(kx) sqrt(kx.^2 + k0(i)^2*e3);
    
    %
    Z1 = @(kx) kz1(kx)./(w(i)*e1*e0);
    Z2 = @(kx) kz2(kx)./(w(i)*e2*e0);
    Z3 = @(kx) kz3(kx)./(w(i)*e3*e0);
    
    %     sigma = Ns*e^2*tau/ms./(1 + 1i*w(i)*tau);
    
    Zup = @(kx) Z2(kx).*(Z1(kx) + 1i*Z2(kx).*tanh(kz2(kx)*d1))./ ...
        (Z2(kx) + 1i*Z1(kx).*tanh(kz2(kx)*d1));
    %     Zup = @(kx) Z1(kx);
    
    Zdown = @(kx) 1i*Z2(kx).*tan(kz3(kx)*d2);
    %     Zdown = @(kx) Z2(kx); % incase of no back gate
    
    F = @(kx) (1./Zup(kx) + 1./Zdown(kx) + 1./sigma(i))/k0(i);
    
    r = newtzero(F,100);
    F(r);
    fval = F(r);
    s = min((fval));
    D = find(F(r) == s);
    if length(D) > 1
        
        if real(r(D(1)) > 0)
            root = vertcat(root,r(D(1)));
        else
            root = vertcat(root,r(D(2)));
        end
    else
        root = vertcat(root,r(D(1)));
    end
    

%     root = -1i*(r);
%     root = vertcat(root, r((1)));
end
plot(imag(root).*k_0', f)
hold on
plot(f/c, f)
% neff = sqrt(ef - (r/k0).^2)

