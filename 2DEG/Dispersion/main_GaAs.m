% Dispersion Relation for 2DEG GaAs
% clear all;
close all;clf
%
% 
% Properties of Graphene Sample taken from:
% Palik

% Define frequency sweep
w = linspace(1e11,1e14, 5e3);
% w = 2*pi*linspace(1e7,10e12, 1e5);


c = 3e8;
ep0_0 = 12.9;
ep0_inf = 11.0;
lambda_l = 1/(292.77e2);
w_l = c/lambda_l;
lambda_t = 1/(268.7e2);
w_t = c/lambda_t;
lambda_gt = 1/(2.4e2);
g_t = c/lambda_gt;
k_t = w_t/c;
tau = 1/g_t;
% InAs
% c = 3e8;
% ep0_0 = 12.9;
% ep0_inf = 11.8;
% w_t = 4.12e13;
% w_l = 4.58e13;
% k_t = w_t/c;
% w_p = 2.95e13;
% W_p = 2.61e12;
% m = .021*m_e;
% Constant Parameters
h = 6.626e-34/(2*pi);
c = 3e8;
e2 = (2.3068e-28);
q = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

% Calculate Plasma Frequency of GaAs
% N = .9e15*1e4;
Ns1 = 1e11*1e4;
m = .063*m_e;

% W_p = (4*pi*Ns*e2/m/c); % Square of plasma frequency
% w_p = sqrt(W_p)


% ep0 = ep0_inf*(w.^2 - w_l^2)./(w.^2 - w_t^2) - w_p^2./w.^2;
% Get Lorentz-based materials
% ep0 = Lorentz(ep0_0, ep0_inf, w_t, g_t,w);
ep0 = ep_GaAs_new(w);
sigma1 = Ns1*e2*tau/m./(1 + 1i*w*tau);
W_p1 = (4*pi^2*Ns1*e2/m); % Square of plasma frequency
w_p1 = sqrt(W_p1);
W_p1 = 4*pi*Ns1*q/(m*c);
w_p1 = sqrt(W_p1);
% 2DEG MAterial
% ep0 = 1;
temp1 = (w/c).^2.*ep0 + (w.^2/W_p1).^2.*ep0;
% temp = (w/c).^2.*ep0 + w.^4.*ep0.^2*m^2/(4*pi^2*(Ns*e2)^2);
% temp = (w/c).^2.*ep0 - ep0.^2.*w.^2./sigma.^2;
k1 = sqrt(temp1);



e0 = 8.85e-12;
ep2 = ep_AlGaAs_new(w);
kx1 = -1i*(ep0 + ep2).*w./sigma1;
kx2 = -1i*(ep0 + 1).*w./sigma1;



%% Dispersion Relation
% close all;clf;figure(1)

% N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% h1 = plot(w*1e-12,real(ep0), 'linewidth',1.4);
% hold on
% h2 = plot(w*1e-12,real(ep2), 'linewidth',1.4);
% plot(w*1e-12,imag(ep0), 'linewidth',.5, 'LineStyle','--')
% plot(w*1e-12,imag(ep2), 'linewidth',.5, 'LineStyle','--')
% set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
% xlim([0 15])
% ylim([-35 55])
% %
% xlabel('$f (\mathrm{THz})$','interpreter','latex')
% ylabel('$\E (\O)$','interpreter','latex')
% legend([h1 h2], {'$\mathrm{GaAs}$', '$\mathrm{Al_{.54}Ga_{.46}As}$'},...
%     'location','northwest','interpreter','latex','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% % grid on
% box on
% % cleanfigure();
% % matlab2tikz('filename',sprintf('epsilon_gaas.tex'));
% % % 
% % %
% % 
% % 
figure(2)
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(real(k1), w*1e-12, 'linewidth',1.4);
hold on
% h2 = plot(real(k2), w, 'linewidth',1.4);

h2 = plot(abs(imag(k1)), w*1e-12, 'linewidth',.5, 'LineStyle',':');
% plot(abs(imag(k2)), w, 'linewidth',.5, 'LineStyle',':');
% 2

h3 = plot(real(w/c)*sqrt(ep0_0), w*1e-12, 'linewidth',1.0, 'color','black');
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
legend([h1 h2 h3],{'$\Re \mathrm{k_x}$',...
    '$\Im \mathrm{k_x}$', 'Light Line'},...
    'location','southeast','interpreter','latex');
xlim([1e3 1e7])
ylim([1e-2 1e1]*w_t)
ylim([1e11 1e14]*1e-12)
xlabel('$k_x (\mathrm{m^{-1}})$','interpreter','latex')
ylabel('$f  (\mathrm{THz})$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on
% % % cleanfigure();
% % matlab2tikz('filename',sprintf('dispersion_gaas.tex'));
% % 
% % 
% % % close all;clf
figure(3)
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(abs(real(kx1)), w*1e-12, 'linewidth',1.4);
hold on
h2 = plot(abs(imag(kx1)), w*1e-12, 'linewidth',1.0, 'LineStyle','--');

% Second Material

% plot(abs(real(kx2)), w, 'linewidth',1.4);
hold on
% plot(abs(imag(kx2)), w, 'linewidth',1.4, 'LineStyle','--');

h3 = plot((w/c)*sqrt(ep0_0), w*1e-12, 'linewidth',.5, 'color','black');
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% 
% xlim([3.5e-2 1e3]*k_t)
% ylim([w(1) 1e14])
xlim([1.2e3 1e10])
xlim([.33e3 1e10])
ylim([1e11 1e14]*1e-12)

xlabel('$k_x (\mathrm{rad/m})$','interpreter','latex')
ylabel('$f (\mathrm{THz})$','interpreter','latex')
legend([h1 h2 h3],{'$\Re \mathrm{k_x}$',...
    '$\Im \mathrm{k_x}$', 'Light Line'},...
    'location','southeast','interpreter','latex');
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on
% % % cleanfigure();
% % matlab2tikz('filename',sprintf('dispersion_gaas_mult.tex'));
% % 
% % 
% % 
% % % Propagation Length
% % % close all;clf
figure(4)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
plot((w)*1e-12, (1./(imag(2*kx2)))*1e6, 'linewidth',1.4);
% symlog()
hold on
% plot(w, abs(1./imag(kx2)), 'linewidth',1.4, 'LineStyle','--');
hold on
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'lin') % Set axes fonts to Times New Roman
% 
% xlim([3.5e-2 1e3]*k_t)
% ylim([w(1) 1e14])
% xlim([30 500]) % 1 THz to 10 THz 
% xlim([30 300]) % 1 THz to 10 THz 
% xlim([30 500]) % 1 THz to 10 THz 
% xlim([.1 10])
% xlim([40e-3 3.5e-1]*1e3)
% ylim([1e-1 1e1])
% ylim([1.7e11 5e14])


ylabel('$L_p (\mathrm{\mu m})$','interpreter','latex')
xlabel('$f(\mathrm{THz})$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on
% % % cleanfigure();
% matlab2tikz('filename',sprintf('plength_gaas_mult_f.tex'));
% % 
% 
% % % Loss TAngent
% figure(5)
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% % Ltan = imag(ep0).*imag(ep2)./(real(ep0).*real(ep2));
% Ltan = (imag(ep0) + imag(ep2))./(real(ep0) + real(ep2));
% Ltan_db = mag2db(abs(Ltan));
% plot(w*1e-12, Ltan_db, 'linewidth',1.4);
% % 
% set(gcf,'Color','white'); % Set background color to white
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
% % % 
% xlim([0e12 10e12]*1e-12)
% % % ylim([w(1) 1e14])
% % % xlim([2e3 1e10])
% % % ylim([1.7e11 5e14])
% % 
% ylabel('$\tan \delta (\mathrm{dB})$','interpreter','latex')
% xlabel('$frequency(\mathrm{Hz})$','interpreter','latex')
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
% box on
% % cleanfigure();
% matlab2tikz('filename',sprintf('ltan_gaas_mult.tex'));


% dBm = 1./(imag(2*kx1));
% % close all;clf
% figure(6)
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% semilogy(w*1e-12, dBm, 'linewidth',1.4);
% hold on
% % plot(w, abs(imag(2*pi./kx1)), 'linewidth',1.4, 'LineStyle','--');
% 
% % Second Material
% 
% % plot(abs(real(kx2)), w, 'linewidth',1.4);
% hold on
% % plot(abs(imag(kx2)), w, 'linewidth',1.4, 'LineStyle','--');
% 
% % plot((w/c.*sqrt(ep0_0)), w, 'linewidth',1.0, 'color','black');
% set(gcf,'Color','white'); % Set background color to white
% % set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% % 
% xlim([0 15])
% % ylim([1.7e11 5e14])
% 
% ylabel('$\frac{1}{2 \Im k_x}  (\mathrm{dB})$','interpreter','latex')
% xlabel('$f(\mathrm{THz})$','interpreter','latex')
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
% box on
% matlab2tikz('filename',sprintf('pl_gaas_mult.tex'));


% wavelength reduction
% close all;clf
% figure(7)
% N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% lambda = c./w;
% lambda_sw = lambda./sqrt(ep0 + ep2);
% plot(lambda*1e6, real(lambda_sw)*1e6, 'linewidth',1.4);
% hold on
% % plot(w, abs(1./imag(kx2)), 'linewidth',1.4, 'LineStyle','--');
% hold on
% set(gcf,'Color','white'); % Set background color to white
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
% % 
% % xlim([3.5e-2 1e3]*k_t)
% % ylim([w(1) 1e14])
% xlim([30 300]) % 1 THz to 10 THz 
% % xlim([5e-2 3e0])
% % xlim([40e-3 3.5e-1]*1e3)
% % ylim([1e-1 1e1])
% % ylim([1.7e11 5e14])
% 
% ylabel('$\lambda_{sp}(\mathrm{\mu m})$','interpreter','latex')
% xlabel('$\lambda(\mathrm{\mu m})$','interpreter','latex')
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% % grid on
% box on
% cleanfigure();
% matlab2tikz('filename',sprintf('wavelength_gaas_mult.tex'));