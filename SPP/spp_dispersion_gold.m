close all
clf
% Classical Interband Drude-Sommerfeld Model of Gold
c = 3e8;
w =  2*pi*linspace(1e14,1.5e15 ,1e3);
% w =  2*pi*linspace(5e14,6.0e14 ,1e4);
f = w/(2*pi);
l = 3e8./f;
k0 = 2*pi./l;

wp = 45e14;
gamma = 8.35e-16;
ll = 450e-9;
w0 = 2*pi*c/ll;
e_d = 4;
e_g = 1 + wp^2./((w0^2 - w.^2) - 1i*gamma*w);

% Drude Critical points model
% Vial and Laroche
% Au 1.1431 1.3202E16 1.0805E14
% 0.26698 ?1.2371 3.8711E15
% 4.4642E14 3.0834 ?1.0968 4.1684E15 2.3555E15
%
% ep_inf = 12.1431; % Check This value
ep_inf = 1.1431;
% ep_inf = 2.63; % 200-nm-thick single-crystalline gold film
wd = 1.3202e16;
gamma = 1.0805e14;

% Data Palik
load gold_JC.mat
N = n + 1i*k;
epsilon = N.^2;
% w = 2*pi*c./lambda;
% f = c./lambda;
% e_g = epsilon;
% k0 = 2*pi./lambda;
% l = lambda;

A1 = .26698;
phi1 = -1.2371;
W1 = 3.8711e15;
T1 = 4.4642e14;

A2 = 3.0834;
phi2 = -1.0968;
W2 = 4.1684e15;
T2 = 2.3555e15;

A = [A1 A2];
phi = [phi1 phi2];
W = [W1 W2];
T = [T1 T2];

e_g = ep_inf - wd^2./(w.^2 + 1i*gamma*w);
for i = 1 : length(A)
    e_g = e_g + A(i)*W(i)*(exp(1i*phi(i))./(W(i) - w - 1i*T(i)) + ...
        exp(-1i*phi(i))./(W(i) + w + 1i*T(i)));
end

% k0 = 2*pi./l;
k_sp = k0.*sqrt(e_g*e_d./(e_g+e_d));
k_x_real = w/c.*sqrt(real(e_g)*e_d./(real(e_g)+e_d));
k_x_imag = w/c.*imag(e_g)./(2*(real(e_g)).^2).*(real(e_g)*e_d./(real(e_g)+e_d)).^(3/2);
beta = real(k_sp);
alpha = imag(k_sp);




figure(1)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(l*1e9, real(e_g), 'linewidth',1.4);
hold on
h2 = plot(l*1e9,imag(e_g), 'linewidth',1.4);
plot(lambda*1e3, real(epsilon),'o','MarkerSize',5);
plot(lambda*1e3, imag(epsilon),'s','MarkerSize',5);

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
matlab2tikz('filename',sprintf('ep_gold.tex'));



figure(2)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(beta*1e-6, f*1e-12, 'linewidth',1.4);
hold on
h2 = plot(alpha*1e-6,f*1e-12, 'linewidth',1.4, 'LineStyle',':');
h3 = plot(sqrt(e_d)*k0*1e-6,f*1e-12,'k',...
    'LineWidth',1)
% h = plot(abs(k_x_real),f);

% plot(abs(k_x_imag),f);
set(gcf,'Color','white');
% xlim([1e6 4e7]*1e-6)
% ylim([1e14 7e14]*1e-12)
xlim([ 4 40])
ylim([105 870])
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
matlab2tikz('filename',sprintf('disp_gold.tex'));