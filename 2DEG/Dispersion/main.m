% Dispersion Relation for 2D Graphene
% Can be modified for 2DEG materials
% 
% Properties of Graphene Sample taken from:
% Terahertz plasmonics in ferroelectric-gated graphene

% Define frequency sweep
w = 2*pi*linspace(0,14e12, 1000);


ep0_0 = 41.5;
ep0_inf = 19.5;
epe_inf = 10;
epe_0 = 26;
w_t0 = 2*pi*4.6e12;
w_te = 2*pi*7.6e12;
g_0 = .51e12;
g_e = .84e12;

ep0_0 = 41.5;
ep0_inf = 19.5;
epe_inf = 10;
epe_0 = 26;
w_t0 = 2*pi*4.6e12;
w_te = 2*pi*7.6e12;
g_0 = .51e12;
g_e = .84e12;

% Constant Parameters
h = 6.626e-34/(2*pi);
c = 3e8;
e = 1.60218e-19; 

% Fermi Levels
E_f = .05 * 1.60218e-19; 
v_f = 1e8*1e-2;

mu = 8e4*1e-4;


% Get Lorentz-based materials
ep0 = Lorentz(ep0_0, ep0_inf, w_t0, g_0,w);
epe = Lorentz(epe_0, epe_inf, w_te, g_e,w);


% 2DEG MAterial



%% Dispersion Relation
close all;clf;figure(1)
temp = epe.*w.^2/c + ep0.*epe/(2*e^2*E_f/h^2)^2.*w.^2.*(w + 1i*e*v_f^2/(mu*E_f)).^2;
temp = ep0.*w.^2/c + ep0/(2*e^2*E_f/h^2)^2.*w.^2.*(w + 1i*e*v_f^2/(mu*E_f)).^2;
k = sqrt(temp);

plot(real(k), w);
hold on
k_c = w/c;
% plot(.5e13*k_c,w)
ylim([0 4e13])