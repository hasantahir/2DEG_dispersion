% function cond = cond_gaas(w)

% From Burke's paper
% High frequency conductivity of the high-mobility two-dimensional
% electron gas
%% Parameters
close all ; clc ; clear all
clf
% w = 2*pi*27.8e12;
% w =  2*pi*linspace(6.2e12,6.4e12, 1e3);
w =  2*pi*linspace(10e12,30e12, 2e3);
q = 1.6012e-19; % Free elecrton charge
m = 9.109e-31; % Free electron mass
h = 6.626e-34;
%% GaAs/AlGaAs
% N = 2e11*1e4; % Electon Density % Sample D
% N = 1.6e11*1e4; % Electon Density % Sample C
% % mu = 15.6*1e6*1e-4; % Mobility % Sample D 
% mu = .6*1e6*1e-4; % Mobility % Sample C
% ms = .063*m; % Effective mass
% % tau = mu*ms/q; % Average Scattering time
% tau = 1e-11;
% ep = 11;
%% GaN/AlGaN From Muravjov's paper
N = 7.5e12*1e4;

% 3 K
mu = 1e6*1e-4;
tau = 1.14e-10;
% 77 K
% mu = 1e4*1e-4;
% tau = 1.14e-12;
% 295 K
% mu = 1200*1e-4;
% tau = .14e-12;
%
ms = .2*m;
ep = 9.5;

% tr = 1e-11;
% Conductivity
cond = N*q^2*tau/ms./(1 - 1i*w*tau);
cond0 = 1;%pi*q^2/(2*h);
% Generate a material for FDTD based simulations
% ep = 12.9;
e_0 = 8.85e-12;
delta = 5e-9;

epsilon =  (1 + 1i*cond./(w*delta*ep*e_0));
% Cond = cond./(epsilon);

% figure(1)
% semilogy(w/(2*pi),real(cond)/cond0,'linewidth',1.4)
% hold on
% semilogy(w/(2*pi),imag(cond)/cond0)


close all;clf;figure(1)
% 
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(w/(2*pi)*1e-12,real(cond), 'linewidth',1.4);
hold on
h2 = plot(w/(2*pi)*1e-12,imag(cond), 'linewidth',1.4);
% plot(w*1e-12,imag(ep0), 'linewidth',.5, 'LineStyle','--')
% plot(w*1e-12,imag(ep2), 'linewidth',.5, 'LineStyle','--')
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% xlim([1 10])
% ylim([-35 55])
%
xlabel('$f (\mathrm{THz})$','interpreter','latex')
ylabel('$\sigma_s (\omega)$','interpreter','latex')
legend([h1 h2], {'Real', 'Imag'},...
    'location','northeast','interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
% % cleanfigure();
% matlab2tikz('filename',sprintf('cond_2DEG_gas_hiT.tex'));

%%%%%%%%%%%
%%%%%%%%%%
 
figure(2)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(w/(2*pi)*1e-12,real(epsilon), 'linewidth',1.4);
hold on
h2 = plot(w/(2*pi)*1e-12,1e2*imag(epsilon), 'linewidth',1.4);
% plot(w*1e-12,imag(ep0), 'linewidth',.5, 'LineStyle','--')
% plot(w*1e-12,imag(ep2), 'linewidth',.5, 'LineStyle','--')
set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',15);
% set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'log') % Set axes fonts to Times New Roman
% xlim([0 20])
% ylim([-30 10])
%
xlabel('$f (\mathrm{THz})$','interpreter','latex')
ylabel('$\varepsilon (\omega)$','interpreter','latex')
legend([h1 h2], {'real', 'imag'},...
    'location','southeast','interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
% % cleanfigure();
% matlab2tikz('filename',sprintf('eps_2DEG_gas.tex'));
% semilogx(w/(2*pi),real(epsilon),'o')
% hold on
% semilogx(w/(2*pi),imag(epsilon),'s')