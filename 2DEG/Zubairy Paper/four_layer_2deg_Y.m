%% solution for dispersion relation for 4 layer 2deg structure

clear all

f = linspace(.1e12, 35e12, 3e1); % Frequency in Hz
% f = 25e12;
% f = 8.3e12; % Frequency
% f = 25e12;
%% Constants
e0 = 8.85e-12; % Free-space permittivity
u0 = 1.25e-06; % free-space permeability
c = 3e8; % speed of light
e = 1.60218e-19; % electron charge
me = 9.109e-31; % Electron mass

% free space wavenumber
w = 2*pi*f;
k0 = 2*pi*f/c;
f0 = 1e12;
k_0 = 2*pi*f0/c;

%% GaAs/AlGaAs herterostructure Burke
% layer material properties 
% Ns =  2e11*1e4; % Electon Density % Sample D
% ms = .063*me; % Effective mass
% % mu = 15.6*1e6*1e-4; % Mobility % Sample D 
% tau = mu*ms/q; % Average Scattering time
% e1 = 1; % Air
% e2 = 11.88; % AlGaAs
% e3 = 13.01; % GaAs

%% GaN/AlGaN herterostructure Moravjov 
% layer material properties
Ns = 5e11 *1e4;
ms = .2*me; % Effective mass
mu = 1e6*1e-4;
tau = 1.14e-10;
e1 = 1; % Air
e2 = 9.5; % AlGaN
e3 = 9.5; % GaN


% Layer geometry
d1 = 20e-9;
d2 = 50e-9;



root = [];
KZ = [];
sigma = Ns*e^2*tau/ms./(1 + 1i*w*tau);
for i = 1 : length(w)
    
    % Ensure decay of the wave at +infinity, we need a negative imaginary
    % part of the kz1
%     kz1 = @(kx) sqrtbr(k0(i)^2 - kx.^2, pi/2);
    kz1 = @(kx) sqrt(k0(i)^2 - kx.^2);
    kz2 = @(kx) sqrt(e2*k0(i)^2 - kx.^2);
    % To change the form of Ydown from cot -> coth
      kz3 = @(kx) sqrt(e3*k0(i)^2 - kx.^2);
%     kz3 = @(kx) 1i*sqrt(-e3*k0(i)^2 + kx.^2);
    
    %
    Z1 = @(kx) kz1(kx)./(w(i)*e1*e0);
    Z2 = @(kx) kz2(kx)./(w(i)*e2*e0);
    Z3 = @(kx) kz3(kx)./(w(i)*e3*e0);
    
    Y1 = @(kx) 1./Z1(kx);
    Y2 = @(kx) 1./Z2(kx);
    Y3 = @(kx) 1./Z3(kx);
    
    Gamma = @(kx) (Z1(kx) - Z2(kx))./(Z1(kx) + Z2(kx)).*exp(-2i*kz2(kx)*d1);
    %     sigma = Ns*e^2*tau/ms./(1 + 1i*w(i)*tau);
    
    Yup = @(kx) Y2(kx).*(1 - Gamma(kx))./(1 + Gamma(kx));
    %     Zup = @(kx) Z1(kx);
    
%     Ydown = @(kx) -Y1(kx).*coth(kz3(kx)*d2);
    Ydown = @(kx) -1i*Y1(kx).*cot(kz3(kx)*d2);
    %     Zdown = @(kx) Z2(kx); % incase of no back gate
    
    F = @(kx) (Yup(kx) + Ydown(kx) + sigma(i));
    
    r = newtzero(F,k0(i));
    % Choose the best value based on absolute tolerance
%     find((abs(F(r))) == min(abs(F(r))))
%     fmin = min(abs(F(r)));
%     fval = abs(F(r));
%     s = min((fval));
% %     root = r(1);
%     D = find(fval == fmin);
%     if length(D) > 1
%         
%         if real(r(D(1)) > 0)
%             root = vertcat(root,r(D(1)));
%         else
%             root = vertcat(root,r(D(2)));
%         end
%     else
%         root = vertcat(root,r(D(1)));
%     end
    
    % choose the maximum real part one
    maxr = find(real(r) == max(real(r)));
    root = vertcat(root, conj(r(maxr)));

%     root = -1i*(r);
%     root = vertcat(root, r((1)));
    KZ = vertcat(KZ, kz1(r(1)));
end
% r(D)
% close all
% plot(abs(real(root)), f)
% hold on
% % plot(abs(imag(root)), f)
% plot(f/c*sqrt(e3), f)
% set(gca,'Xscale','Lin','Yscale','Lin')
% xlim([1e4 1e6])


close all;clf;figure(1)
% 
N = 3;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(abs(real(root)), f*1e-12, 'linewidth',1.4);
hold on
% h3 = plot(abs(imag(root)), f*1e-12, 'linewidth',1.4);
h2 = plot(f/c*10, f*1e-12, 'linewidth',1.0);
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
xlim([3e4 1.1e7])
ylim([.1 35])
%
xlabel('$k_x (\mathrm{rad~m}^{-1})$','interpreter','latex')
ylabel('$f (\mathrm{THz})$','interpreter','latex')
legend([h1 h2], {'Back gated', 'Light Line x 10'},...
    'location','southeast','interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
% % % cleanfigure();
% matlab2tikz('filename',sprintf('gated_disp.tex'));
