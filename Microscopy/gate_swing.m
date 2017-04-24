% This code shows tunability of plasma resonance with gate bias
Lg1 = .05e-6;
Lg2 = .01e-6;
q = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass
e0 = 8.85e-12;

ep = 12.7;
ms = .067*m_e;
dg = 20e-9;
N0 = 2e12*1e4;

Vt = -0.764;
V = linspace(0 , 5, 1e3);

Ns = N0 * ( 1 - V/Vt);

omega1 = pi/(2*Lg1)*sqrt(q^2*dg*Ns/(e0*ep*ms));
omega2 = pi/(2*Lg2)*sqrt(q^2*dg*Ns/(e0*ep*ms));


plot(V,omega1)
% shading interp
close all;clf;figure(1)
% 
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(V,omega1*1e-12, 'linewidth',1.4);
hold on
h2 = plot(V,omega2*1e-12, 'linewidth',1.4);
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'log','XScale', 'lin') % Set axes fonts to Times New Roman

%
xlabel('$V_g (\mathrm V)$','interpreter','latex')
ylabel('$f (\mathrm{THz})$','interpreter','latex')

legend([h1 h2], {'$L_g = 1 \mu m$', '$L_g = 100 nm$'},...
    'location','northwest','interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on

%% figure(2)
% shading interp
figure(2)
% 
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(V,Ns*1e-4, 'linewidth',1.4);
set(gcf,'Color','white');
set(gca,'FontName','times new roman','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman

%
xlabel('$V_g (\mathrm V)$','interpreter','latex')
ylabel('$N_s (\mathrm{cm}^{-2})$','interpreter','latex')

% legend([h1 h2], {'$L_g = 1 \mu m$', '$L_g = 100 nm$'},...
%     'location','northwest','interpreter','latex','FontSize',15);
set(gca,'FontName','times new roman','FontSize',15) % Set axes fonts to Times New Roman
% grid on
box on
