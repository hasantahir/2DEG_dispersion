function [sig] = Sigma_Graphene( g, G, M, T, lam )
% Surface conductivity of graphene (RPA model)
% IN:   g   - scattering rate for intraband part, [eV], e.g. G = 11e-3
%       G   - scattering rate for interband part, [eV], e.g. G = 11e-3
%       M   - chemical potential, [eV], typically m = 0.1 : 0.6
%       T   - temperature, [eV], e.g. T=300*kB (room temperature)
%       lam - range of wavelength, [m]
% OUT:  sig - complex array, surface conductivity of graphene

G0 = 1e-5;
if abs(G)<G0, fprintf('PV-version is used\n'); end

% Constants ---------------------------------------------------------------
c  = 299792458;                 % speed of light, [m/s]
kB = 8.61733247800e-05;         % Boltzmann constant, [eV/K]
h_ = 6.58211928150e-16;         % reduced Planck constant, [eV*s]
e  = 1.60217656535e-19;         % electron charge, [C]
s0 = e/4/h_;                    % conductivity of T/m=+0,m=+0 graphene, [S]
e0 = 8.854187817620e-12;        % permittivity in vacuum, [F/m] 

w = h_*2*pi*c./lam;                               % frequencies, [eV]
f = @(a) (tanh(0.5*(a-M)/T) + tanh(0.5*(a+M)/T)); % integral kernel

% Carrier density to chemical potential conversion
% vf = 1e6; 						% Fermi velocity of electrons, [m/s]
% ns = 1e13*1e4; 					% sample carrier density, [1/m^2]
% m = h_*vf*sqrt(pi*ns); 			% chemical potential, derived from carrier density (if ns=1e131/cm^2 then m ~ 0.36893 eV)

% Intraband Conductivity, [S] ---------------------------------------------
sig1 = 8i*s0/pi * log(2*cosh(0.5*M/T)) * T./(w+1i*g);

% Interband Conductivity, [S] ---------------------------------------------
if abs(G)<G0, func = @(W) quadgk( @(a) ( f(a)-f(W/2) )./(4*a.^2-W.^2       ), 0,inf, 'AbsTol',0);
else          func = @(W) quadgk( @(a)   f(a)         ./(4*a.^2-(W+1i*G).^2), 0,inf, 'AbsTol',0);    
end
integral = arrayfun(func,w);
sig2 = - 2i*s0/pi*(w+1i*G).*integral;

% if G~0 we have real PV integral => needs to add imaginary part
if abs(G)<G0, sig2 = sig2 + s0/2*f(w/2); end

sig = sig1 + sig2;

end
