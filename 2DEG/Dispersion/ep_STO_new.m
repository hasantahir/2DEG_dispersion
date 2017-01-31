function [epr,epi] =  ep_STO_new(w)
% Create Dielectric Constant of SrTi2O3
% Based on M. Boucherit
% Extreme charge density SrTiO3/GdTiO3 heterostructure 
% field effect transistors
lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
ep1 = 1;

%% Data from Paper
% N = 1e14*1e4; % Carrier density in cm^-2
% mu = 6*1e-4; % Mobility cm^2 /V/s

% w = linspace(1e11,1e14, 1e5);
ep_inf = 5.2;


% N = 3e14*1e4;
% mu = 100*1e-4;

% N = 6e13*1e4;
% mu = 100*1e-4;

% m = 9.109e-31; % Electron mass
% ms = 6*m; % Effective mass of electron
% e2 = (2.3068e-28);

% Relaxation time
l1 = 18.4e-6;
w1 = 544e2;
g1 = .049*w1;
cons1 = 1.56;

l2 = 56.3e-6;
w2 = 178e2;
g2 = .039*w2;
cons2 = 3.6;

l3 = 114.3e-6;
w3 = 87.7e2;
g3 = .5*w3;
cons3 = 311;

wi = c*[ w1 w2 w3];
gi = c*[ g1 g2 g3];
consi = [ cons1 cons2 cons3];

epr = zeros(size(w));
epi = zeros(size(w));

epr = epr + ep_inf;
epi = 0;

for i = 1 : length(wi)
    epr = epr + consi(i)*wi(i)^2*(wi(i)^2 - w.^2)./...
        ((wi(i)^2 - w.^2).^2 + gi(i)^2*w.^2);
    
    epi = epi + consi(i)*wi(i)^2*(gi(i)*w)./...
        ((wi(i)^2 - w.^2).^2 + gi(i)^2*w.^2);
end

% close all;clf;figure(1)
% N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% hold on
% % plot(w,epr, 'linewidth',1.4)
% % close all;clf;figure(2)
% % plot(w,epi, 'linewidth',1.4, 'LineStyle','--')
% plot(c./(w)*1e6, epi, 'linewidth',1.4);
% set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% % xlim([6e12 10e12])
% xlim([0 500])
% % ylim([1e-2 1e2])
% %
% xlabel('$f (\mathrm{Hz})$','interpreter','latex')
% legend({'$\mathrm{\varepsilon_{real}}$', '$\mathrm{\varepsilon_{imag}}$'},...
%     'location','northwest','interpreter','latex','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% % grid on
% box on       
    
    
    
    
    
  