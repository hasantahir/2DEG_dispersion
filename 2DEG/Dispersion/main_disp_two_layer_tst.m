close all;clear all
% Define frequency sweep
w = 2*pi*linspace(1e8,1e15, 1e2);
% w = 2*pi*linspace(1e7,10e12, 1e5);


c = 3e8;
ep0_0 = 12.9;
ep0_inf = 11.0;
lambda_l = 1/(292.77e2);
w_l = 2*pi*c/lambda_l;
lambda_t = 1/(268.7e2);
w_t = 2*pi*c/lambda_t;
lambda_gt = 1/(2.4e2);
g_t = 2*pi*c/lambda_gt;
k_t = w_t/c;

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
% e = 1.60218e-19; 
m_e = 9.109e-31; % Electron mass

% Calculate Plasma Frequency of GaAs
% N = .9e15*1e4;
Ns = 5e12*1e4;
m = .063*m_e;

% W_p = (4*pi*Ns*e2/m/c); % Square of plasma frequency
% w_p = sqrt(W_p)

root = [];
% ep0 = ep0_inf*(w.^2 - w_l^2)./(w.^2 - w_t^2) - w_p^2./w.^2;
% Get Lorentz-based materials
for i = 1 : length(w)
ep2 = Lorentz(ep0_0, ep0_inf, w_t, g_t,w(i));
sigma = Ns*e2./(1i*m*w(i));
W_p = (4*pi*Ns*e2/m); % Square of plasma frequency
w_p = sqrt(W_p);


ep1 = 1;
k1 = (w(i)/c);
k2 = (w(i)/c).*sqrt(ep2);

kz1 = @(kx) sqrt(k1.^2 - kx.^2);
kz2 = @(kx) sqrt(k2.^2 - kx.^2);
func = @(kx) (ep2./kz1(kx) + 1./kz2(kx) - 1i*Ns*e2./(m.*w(i).^2));


% p = linspace(lxlim,uxlim,num);
% root = [];
% for i = 1 : length(p)
    r = newtzero(func,k2);
    root = vertcat(root,r);
end
% figure(2)
N = 2;
axes('ColorOrder',brewermap(N,'Set1'),'NextPlot','replacechildren')
h1 = plot(real(root(1:length(w))), w/(2*pi), 'linewidth',1.4);
hold on
h2 = plot(imag(root(1:length(w))), w/(2*pi), 'linewidth',1.4);

set(gcf,'Color','white'); % Set background color to white
% set(gca,'FontName','times new roman','FontSize',11,'YScale', 'log','XScale', 'log') % Set axes fonts to Times New Roman
% legend([h1 h2],{'$N_s = 5 \times 10^{12} \mathrm{cm^{-2}}$',...
%     '$N_s = 5 \times 10^{15} \mathrm{cm^{-2}}$'},...
% %     'location','southeast','interpreter','latex');
% xlim([3.5e-2 1e3]*k_t)
% ylim([1e-2 1e1]*w_t)
% ylim([1e-2 1e1]*w_t)
xlabel('$k (\mathrm{rad/m})$','interpreter','latex')
ylabel('$\omega (\mathrm{rad/s})$','interpreter','latex')
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
grid on
box on

