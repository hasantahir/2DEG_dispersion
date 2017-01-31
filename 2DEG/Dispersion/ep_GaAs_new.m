function epr =  ep_GaAs_new(w)
% Palik Vol. 1 pg 431

lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
ep1 = 1;

%% Data from Paper
% N = 1e14*1e4; % Carrier density in cm^-2
% mu = 6*1e-4; % Mobility cm^2 /V/s

% w = linspace(1e12,1e14, 1e5);
ep_inf = 11.0;


wt1 = 268.8e2;
wl1 = 292.1e2;
g1 = 2.4e2;



wt = c*wt1;
wl = c*wl1;
g = c*g1;


epr = zeros(size(w));
% epi = zeros(size(wt));

epr = epr + ep_inf;
epi = 0;

for i = 1 : length(wt)
    epr = epr.*(1+ (wl(i)^2 - wt(i).^2)./...
        (wt(i)^2 - w.^2 - 1i*w*g(i)));
end

% close all;clf;figure(1)
% N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% hold on
% plot(c./(w)*1e6,real(epr), 'linewidth',1.4)
% % close all;clf;figure(2)
% plot(c./(w)*1e6,imag(epr), 'linewidth',1.4, 'LineStyle','--')
% set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'log') % Set axes fonts to Times New Roman
% % xlim([6e12 10e12])
% % ylim([1e-2 1e2])
% xlim([0 100])
% %
% xlabel('$f (\mathrm{Hz})$','interpreter','latex')
% legend({'$\mathrm{\varepsilon_{real}}$', '$\mathrm{\varepsilon_{imag}}$'},...
%     'location','northwest','interpreter','latex','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% % grid on
% box on       
    
    
    
    
    
  