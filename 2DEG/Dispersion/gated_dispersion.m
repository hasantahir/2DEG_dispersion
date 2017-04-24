%% Gated dispersion constant
e_0 = 8.85e-12;
h = 6.626e-34/(2*pi);
c = 3e8;
e2 = (2.3068e-28);
q = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass


Ns = 2e11 *1e4;

ep1 = 9.6;
ep2 = 11;
ms = .063*m_e; % Effective mass

d_in = 20e-9;
d_sc = 50e-9;
% Define a vector for the wavenumber
w = linspace(1e2,1e8, 2e3);

omega = 4*pi*Ns*e^2*q/ms./(ep2*coth (q*d_sc) + ep1*coth(q*d_in));
omega = sqrt(omega)

plot(q,omega)