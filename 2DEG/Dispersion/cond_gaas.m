function cond = cond_gaas(w)

% From Burke's paper
% High frequency conductivity of the high-mobility two-dimensional
% electron gas
%% Parameters
% close all ; clc ; clf
% w = 2*pi*linspace(1e12,10e12, 1e3);
q = 1.6012e-19; % Free elecrton charge
m = 9.109e-31; % Free electron mass

N = 2e11*1e4; % Electon Density % Sample D
% N = 1.6e11*1e4; % Electon Density % Sample C
mu = 15.6*1e6*1e-4; % Mobility % Sample D 
% mu = .6*1e6*1e-4; % Mobility % Sample C
ms = .063*m; % Effective mass
tr = mu*ms/q; % Average Scattering time
% tr = 1e-12;
% Conductivity
cond = N*q^2*tr/ms./(1 - 1i*w*tr);

e_0 = 8.85e-12;
delta = 5e-9;
ep0_0 = 12.9;
ep_inf = 11.0;

epsilon =  (1 + 1i*cond./(w*delta*ep_inf*e_0));
% cond = cond./(epsilon);

% sigma
% cond = 1i*N*q^2./(ms*w);

% Generate a material for FDTD based simulations
% ep = 12.9;
% e_0 = 8.85e-12;
% delta = 3e-9;
% 
% epsilon =  ep + 1i*cond./(w*delta*e_0);

% semilogy(w/(2*pi),real(cond))
% hold on
% semilogy(w/(2*pi),imag(cond))
% 
% figure(2)
% semilogx(w/(2*pi),real(epsilon))
% hold on
% semilogx(w/(2*pi),imag(epsilon))