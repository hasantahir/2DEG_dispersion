function [epr,epi] =  ep_LAO_new(w)
% Microwave Dielectric Properties of Lanthanum
% Aluminate Ceramics and Single Crystal
% Takeshi Shimada, Ken-ich Kakimoto* and Hitoshi Ohsato*
% R&D Center, NEOMAX Co., LTD.
lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
ep1 = 1;

%% Data from Paper
% N = 1e14*1e4; % Carrier density in cm^-2
% mu = 6*1e-4; % Mobility cm^2 /V/s

% w = linspace(1e12,1e14, 1e5);
ep_0 = 310;
ep_inf = 4.0;


% N = 3e14*1e4;
% mu = 100*1e-4;

% N = 6e13*1e4;
% mu = 100*1e-4;

% m = 9.109e-31; % Electron mass
% ms = 6*m; % Effective mass of electron
% e2 = (2.3068e-28);

% Relaxation time
wt1 = 182.8e2;
wl1 = 276.2e2;
gt1 = 5.4e2;
gl1 = 4.1e2;

wt2 = 426.5e2;
wl2 = 428.5e2;
gt2 = 5.4e2;
gl2 = 114.9e2;

wt3 = 430.4e2;
wl3 = 596.3e2;
gt3 = 131.8e2;
gl3 = 12.8e2;

wt4 = 647.3e2;
wl4 = 743.8e2;
gt4 = 39.4e2;
gl4 = 11.6e2;

wt = c*[ wt1 wt2 wt3 wt4];
wl = c*[ wl1 wl2 wl3 wl4];
gt = c*[ gt1 gt2 gt3 gt4];
gl = c*[ gl1 gl2 gl3 gl4];


epr = zeros(size(w));
% epi = zeros(size(wt));

epr = epr + ep_inf;
epi = 0;

for i = 1 : length(wt)
    epr = epr.*(wl(i)^2 - w.^2 - 1i*w*gl(i))./...
        (wt(i)^2 - w.^2 - 1i*w*gt(i));
end
% for i = 1 : length(wi)
%     epr = epr + Si(i)*wi(i)^2*(wi(i)^2 - w.^2)./...
%         ((wi(i)^2 - w.^2).^2 + gi(i)^2*w.^2);
%     
%     epi = epi + Si(i)*wi(i)^2*(gi(i)*w)./...
%         ((wi(i)^2 - w.^2).^2 + gi(i)^2*w.^2);
%     
% end
% Eg =  5.6*1.60218e-19;% eV;
% A = 243.7;
% E0 = 9.025*1.60218e-19;% eV;;
% Gam = 3.989*1.60218e-19;% eV;;
% h = 6.62607004e-34;
% E = h*w;
% epi = A*E0*Gam*(E - Eg).^2./(E.*((E.^2 - Eg^2).^2 + Gam^2*E.^2));

% close all;clf;figure(1)
% N = 2;
% axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
% hold on
% plot(w,real(epr), 'linewidth',1.4)
% % close all;clf;figure(2)
% plot(w,imag(epr), 'linewidth',1.4, 'LineStyle','--')
% set(gcf,'Color','white');
% set(gca,'FontName','times new roman','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'log') % Set axes fonts to Times New Roman
% % xlim([6e12 10e12])
% % ylim([1e-2 1e2])
% %
% xlabel('$f (\mathrm{Hz})$','interpreter','latex')
% legend({'$\mathrm{\varepsilon_{real}}$', '$\mathrm{\varepsilon_{imag}}$'},...
%     'location','northwest','interpreter','latex','FontSize',11);
% set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% % grid on
% box on       
    
    
    
    
    
  