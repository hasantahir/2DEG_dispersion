close all
clf
% Classical Interband Drude-Sommerfeld Model of Gold
c = 3e8;
w =  2*pi*linspace(1e14,1.5e15 ,1e3);
f = w/(2*pi);
l = 3e8./f;
% wp = 45e14;
% gamma = 8.35e-16;
% ll = 450e-9;
% w0 = 2*pi*c/ll;
e_d = 4;
% e_g = 1 + wp^2./((w0^2 - w.^2) - 1i*gamma*w); 
% Ag 15.833 1.3861E16 4.5841E13 
% 1.0171 ?0.93935 6.6327E15 
% 1.6666E15 15.797 1.8087 9.2726E17 2.3716E17
% Drude Critical points model
% Vial and Laroche
% 
ep_inf = 15.833;
% ep_inf = 1.4447;

wd = 1.3861e16;
gamma = 4.5841e13;

load silver_JC.mat
N = n + 1i*k;
epsilon = N.^2;


A1 = 1.1071;
phi1 = -.93935;
W1 = 6.6327e15;
T1 = 1.6666e15;

A2 = 15.797;
phi2 = 1.8087;
W2 = 9.2726e17;
T2 = 2.3716e17;

A = [A1 A2];
phi = [phi1 phi2];
W = [W1 W2];
T = [T1 T2];

e_s = ep_inf - wd^2./(w.^2 + 1i*gamma*w);
for i = 1 : length(A)
    e_s = e_s + A(i)*W(i)*(exp(1i*phi(i))./(W(i) - w - 1i*T(i)) + ...
        exp(-1i*phi(i))./(W(i) + w + 1i*T(i)));
end

k0 = 2*pi./l;
k_sp = k0.*sqrt(e_s*e_d./(e_s+e_d));
k_x_real = k0.*sqrt(real(e_s)*e_d./(real(e_s)+e_d));
k_x_imag = k0.*imag(e_s)./(2*real(e_s).^2).*sqrt(real(e_s)*e_d./(real(e_s)+e_d));
beta = real(k_sp);
alpha = imag(k_sp);




figure(1)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(l*1e9, real(e_s), 'linewidth',1.4);
hold on
h2 = plot(l*1e9,imag(e_s), 'linewidth',1.4)save ;
plot(lambda*1e9, real(epsilon),'o','MarkerSize',5);
plot(lambda*1e9, imag(epsilon),'s','MarkerSize',5);

set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11,'YScale', 'lin','XScale', 'lin') % Set axes fonts to Times New Roman
legend([h1 h2],{'$\Re \mathrm{\E(\O)}$',...
    '$\Im \mathrm{\E(\O)}$'},...
    'location','southeast','interpreter','latex');
xlabel('$\lambda (\mathrm{nm})$','interpreter','latex')
ylabel('$\E(\O)$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% grid on
box on
xlim([188 1000])
% ylim([1e14 7e14]*1e-12)
% % % cleanfigure();
matlab2tikz('filename',sprintf('ep_silver.tex'));



figure(2)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(beta*1e-6, f*1e-12, 'linewidth',1.4);
hold on
h2 = plot(alpha*1e-6,f*1e-12, 'linewidth',1.4, 'LineStyle',':');
h3 = plot(sqrt(e_d)*k0*1e-6,f*1e-12,'k',...
    'LineWidth',1);
% h = plot(real(k_x_real),f);

% plot(imag(k_x_imag),f);
set(gcf,'Color','white');
xlim([2 150])
ylim([1e14 1.2e15]*1e-12)
% set(h,'Color','black','LineWidth',1.4)
% set(gca,'FontSize',10, 'FontName' , 'times new roman','YScale', 'log','XScale', 'log')
xlabel('$k_x (\mathrm{\mu m^{-1}})$',...
    'Interpreter','latex');

ylabel('$f (\mathrm{THz})$',...
    'Interpreter','latex');

legend([h1 h2],{'$\Re \mathrm{k_x}$',...
    '$\Im \mathrm{k_x}$'},...
    'location','southeast','interpreter','latex');
set(gca,'box','on')
hold off
matlab2tikz('filename',sprintf('disp_silver.tex'));