% Dispersion from Nakayama's paper
% InAs
close all;clc;clear all

% Constant Parameters
h = 6.626e-34/(2*pi);
c = 3e8;
e2 = (2.3068e-28);
q = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass


ep_inf = 11.8;
ep_0 = 15;
w_t = 4.18e13;
w_l = 4.58e13;
w_p = 2.95e13;
N = 7e16*1e6;
Ns = 5e13;
m =.021*m_e;

W = 2.61e14;

w = linspace(1e-2, 1e2, 2e4)*w_t;

% Drude model for the dielectric of InAs
ep = ep_inf * (w.^2 - w_l^2)./(w.^2 - w_t^2) - w_p^2./w.^2;

temp = (w/c).^2.*ep.*(1 + (2*w/W).^2.*ep);

k = sqrt(temp);
k_t = w_t/c;

%% 
n = 3;
figure(1)
axes('ColorOrder',brewermap(n,'Set1'),'NextPlot','replacechildren')
h1 = plot(real(k), w, 'linewidth',1.4);
hold on
h2 = plot(imag(k), w, 'linewidth',1.4,'LineStyle','--');
h3 = plot((w/c)*sqrt(ep_inf), w, 'linewidth',1.0,'LineStyle',':');
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% xlim([6e3 1e7])
% ylim([1e-2 1e1]*w_t)
% ylim([1e11 1e13]*1e-12)
xlabel('$k / k_t$','interpreter','latex')
ylabel('$\omega / omega_t $','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on