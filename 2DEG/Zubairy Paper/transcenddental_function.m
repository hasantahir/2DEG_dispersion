w = 2*pi*25e12;
e_0 = 8.85e-12;
c = 3e8;
e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

k0 = w/c;
% Material properties
Ns = 4.14e11 *1e4;
e_d = 12.8;
ms = .067*m_e; % Effective mass
Nmax = 1e5 ; %Number of iterations

kz = linspace(2e5,7e7,1e4);
% layer thicknesses
d = 20e-9;

kxd = sqrt(k0^2*e_d - kz.^2);
kxa = sqrt(k0^2 - kz.^2);

k0 = 2*pi/la0; 
pc = ef/ec; 
ps = ef/es;

s = 1 - 2*mode;      % s = 1,-1 for TM0,TM1

kzm = sqrt(be0.^2 - ef);
kza = sqrt(be0.^2 - ec);
kzs = sqrt(be0.^2 - es);
B = (pc*kza - ps*kzs)/2; 

for N=1:Nmax
   cth = coth(2*k0*kzm.*a);
   A = -kzm.*cth + s * sqrt(B.^2 + kzm.^2.*(cth.^2-1));
   B = sqrt(kzm.^2 + 2*kzm.*A.*cth + A.^2);
   kza = (A+B)/pc;
   ga_new = sqrt(kza.^2 + ec-ef);
   kza = sqrt(ga_new.^2 + ef-ec);            % redundant
   kzs = sqrt(ga_new.^2 + ef-es);
   B = (pc*kza-ps*kzs)/2;
   if norm(ga_new-kzm)<tol, break; end       % break out before N=Nmax
   kzm = ga_new;
end
% kz = linspace(
Gamma = -1i*e*Ns*kxd/(2*ms*e_0*e_d)/w^2;
zeta = kxd./(kxa*e_d);


lw = 'linewidth';
f1 = (1 - Gamma).*(1 + zeta).*exp(-1i*kxd*d);
 
f2 = Gamma.*(1 - zeta).*exp(1i*kxd*d);


close all;clf
plot(kz, abs(real(f1) - real(f2)), lw, 1.5); 
hold on
% plot(kz, abs(real(f2)), lw, 1.5);
set(gca,'Xscale','Log','Yscale','Log')

figure(2)
plot(kz, abs(imag(f1) - imag(f2) ), lw, 1.5); 
hold on
% plot(kz, abs(imag(f2)), lw, 1.5);
set(gca,'Xscale','Log','Yscale','Log')