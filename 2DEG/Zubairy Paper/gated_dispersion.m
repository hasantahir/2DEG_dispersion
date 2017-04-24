% close all; clf;clear
%% Gated dispersion constant
e_0 = 8.85e-12;
h = 6.626e-34/(2*pi);
c = 3e8;
e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass


Ns = 4.14e12 *1e4;

ep1 = 12.8;
ep2 = 3.7;
ms = .067*m_e; % Effective mass

d_in = 90e-6;
% Define a vector for the wavenumber
K = linspace(1/c,50/c, 2e4);

% Ando's expression
omega = Ns*e^2*K/ms./(ep2 + ep1*coth(K*d_in));
omega = sqrt(omega);

%% Chaplik gated expression
%   \O_p^2 = \frac{N_s e^2 k}{m^{\ast} \E}\left(1 + \frac{\E - 1}{\E + 1}\e^{-2k  k d}\right)

% ep2 = 9.5;
% d = 20e-9;
% omega = 2*pi*Ns*e2*q/(ms*ep2).*(1 + (ep2 - 1)/(ep2 + 1).*exp(-2*q*d));
% omega = sqrt(omega);


%% Expression from the Eguiluz
% Ns = 2e12 *1e4;
% ep0 = 3.7;
% eps = 12;
% epm = -1e6;
% ms = .2*m_e;
% d = 1e-5 * 1e-2;
% Bs = sqrt(q.^2 - omega.^2*eps);
% B0 = sqrt(q.^2 - omega.^2*ep0);
% Bm = sqrt(q.^2 - omega.^2*epm);
% omega_p = 8e15;
% gamma = (Bm*ep0 + B0*em).^-1 .*(Bm*ep0 - B0*em);

clf
% omega = z
plot(c*K,omega*c)
hold on
% plot(c*K,K*sqrt(ep1))
% set(gca,'Xscale','Log','Yscale','Log')