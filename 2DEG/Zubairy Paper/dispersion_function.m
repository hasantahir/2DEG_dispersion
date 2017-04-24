e_0 = 8.85e-12;
h = 6.626e-34/(2*pi);
c = 3e8;
e2 = (2.3068e-28);
e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

% Material properties
Ns = 1e12 *1e4;
ep1 = 9.6;
ep2 = 11;
ms = .063*m_e; % Effective mass
mu = 1e6*1e-4;
tau = 1.14e-10;

% layer thicknesses
d0 = 20e-9;
d2 = 50e-9;

% Chinese Homotopy Method
epa = 7.34 ; % GaN/AlGaN layers combined
epb = 11; % Silicon base

% Define a vector for the wavenumber
w = 2*pi*25e12;

% Surface conductivity
sigma = Ns*e^2*tau/ms./(1 - 1i*w*tau);

syms kx
% Input impedance
Z0 = @(kx) 1i*kx./(w*e_0);
Za = @(kx) 1i*kx./(w*epa*e_0);
Zb = @(kx) 1i*kx./(w*epb*e_0);

clf;
% function
f = @(kx) Zb(kx) + 1i*Z0(kx).*tanh(kx*d0) - 1i*coth(kx*d2).*(Z0(kx) + 1i*Zb(kx).*tanh(kx*d0)) + ...
    Zb(kx).*sigma.*(Z0(kx) + 1i*Zb(kx).*tanh(kx*d0));
kx = linspace(0,18,1e3);
plot(kx, real(f(kx)))
hold on;plot(kx, imag(f(kx)))
set(gca,'Xscale','Log','Yscale','Log')
plot(kx,c*kx);



