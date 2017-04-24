clear all;
la0 = 1550e-9;
% ef = 3.5^2;
% es = 1.45^2;
% ef = 3.3^2;
% es = 3.256^2;

%% Kekatpure MDM
% es = -95.92-1i*10.97;
% ef = 2.1025;
% ec = -143.49 - 1i*9.52;
% h = 3e-6;
% h = 300e-9;

%% Kekatpure MDM
% ec = 1;
es = 2.1025;
ef = -143.49 - 1i*9.52;
ec = 2.1025;
h = 100e-9;
% h = 300e-9;

%% Definitions
k0 = 2*pi./la0;
p = ef/ec;
q = ef/es;
Kc = k0*sqrt(ef - ec);
Ks = k0*sqrt(ef - es);
% Qc = k0*sqrt(ec - ef);
% Qs = k0*sqrt(es - ef);

ac = @(k) sqrt(Kc.^2 + k.^2);
as = @(k) sqrt(Ks.^2 + k.^2);
gc = @(k) sqrt(Kc.^2 - k.^2);
gs = @(k) sqrt(Ks.^2 - k.^2);

Gc = @(k) sqrt(k.^2 + p^2*gc(k).^2);
Gs = @(k) sqrt(k.^2 + q^2*gs(k).^2);
S = @(k) (p*ac(k) + q*as(k))/2;

% Kekatpuere defnition (10) and (9)
% f = @(k) tan(k*h/2) - ((p*q*gc(k) .* gs(k) - k.^2) - Gc(k).*Gs(k))./(k.*(p*gc(k) + q*gs(k)));

% Kekatpure MDM

% Plasmonic
% fp = @(k) (k.^2 + 2*S(k).*k.*coth(k*h) + p*q*ac(k).*as(k))/k0;

% Kekatpure DMD
% Plasmonic
fp = @(k) (tanh(k*h) + (k.*(p*ac(k) + q*as(k)))./(k.^2 + p*q*ac(k).*as(k)))/k0;

% Dielectric
f1 = @(k) (tan(k*h) - (k.*(p*gc(k) + q*gs(k)))./(k.^2 - p*q*gc(k).*gs(k)));

% % Orfanidis definition
% la0 = 632.8;
% ef  = (0.0657 - 4i)^2;
% es = 1.5^2;
% ec = 1.55^2;
% p = ef/ec;
% q = ef/es;
% a = 20/2;
% k0 = 2*pi/la0;
% kp = @(beta) sqrt(beta.^2 - ef);
% ac = @(beta) sqrt(beta.^2 - ec);
% as = @(beta) sqrt(beta.^2 - es);
% f2 = @(beta) tanh(2*kp(beta)*a) + (kp(beta).*(p*ac(beta) + q*as(beta)))./(kp(beta).^2 + p*q*ac(beta).*as(beta));
% f2 = @(x) sin(x) + x - 1;
% r = newtzero(f2,1)
% initial_guess = k0;
% p = linspace(lxlim,uxlim,num);
r = newtzero(f1,Ks);
root = r
neff = sqrt(ef - (r/k0).^2)
% kp = kp(r)/k0