close all
clear all
% % Load material data from Gold
load gold_palik.mat
% 
lam = lambda;
n_g = n;
k_g = k;
N = n_g + 1i*k_g;
e_g = N.^2;
c = 3e8;
f = c./lam;
k0 = 2*pi./lam;

wp = 45e14;
gamma = 8.35e-16;
ll = 450e-9;
w0 = 2*pi*c/ll;
e_d = 4;
e_g = 1 + wp^2./((w0^2 - w.^2) - 1i*gamma*w); 
% 
ed = 1;
k_sp = k0.*sqrt(e_g*ed./(e_g+ed));
% 
h = plot(real(k_sp),f);
hold on
plot(imag(k_sp),f);
set(gca,'FontSize',10, 'FontName' , 'times new roman','YScale', 'log','XScale', 'log')
% rng('default'); rng(0)
% data = data + 1e-1*randn(size(data));
% ff = @(x) exp(x).*cos(10*x).*tanh(4*x);
grid = linspace(min(lam),max(lam),length(lam))';
data1 = real(e_g);
data2 = imag(e_g);

%set up for both fit commands in the stats toolbox
% xdata = [-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9];
% ydata = [0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001];

fun = @(p,grid) p(1)*cos(p(2)*grid)+p(2)*sin(p(1)*grid);

pguess = [1,0.2];

lam_cheb = chebfun(grid,[min(lam),max(lam)],'equi','eps',1e-3);
f_nr = chebfun(data1,[min(lam),max(lam)],'eps',1e-4);
f_ni = chebfun(data2,[min(lam),max(lam)],'eps',1e-4);
% fexact = chebfun(real(n_g));
% error = norm(f_nr-fexact,inf)
plot(f_nr)
hold on
plot(f_ni)
plot(lam,real(e_g),'o')

plot(lam,imag(e_g),'o')

N_cheb = f_nr - 1i*f_ni;
e_g_cheb = f_nr + 1i*f_ni;
% k0 = 2*pi./lam_cheb;
% k_sp_cheb = k0.*sqrt(e_g_cheb*e_d./(e_g_cheb+e_d));
% figure(4)
% h = plot(real(k_sp_cheb));
% hold on
% plot(imag(k_sp_cheb));
% % xlim([0 1e8])



figure(2)
x = grid;
y_rng = real(n_g);
y_ing = imag(n_g);
lambda = linspace(min(lam),max(lam),1e1*length(lam));
xq = lambda(2:end-1);
f = c./lambda;
k0 = 2*pi./lambda;
y_nr = spline(x,y_rng,lambda);
y_ni = spline(x,y_ing,lambda);
y_nr_1 = interp1(x,y_rng,lambda,'pchip');
y_ni_1 = interp1(x,y_ing,lambda,'pchip');
plot(x,y_rng,'o',lambda,y_nr_1)
hold on
plot(x,y_rng,'o',lambda,y_ni_1)


% plot(x,y_ing,'o',lambda,y_ni)
% 
% 
figure(3)
N = y_nr - 1i*y_ni;
N = y_nr_1 - 1i*y_ni_1;
eg = (N).^2;
k_sp = k0.*sqrt(eg*ed./(eg+ed));
k_x_real = k0.*sqrt(real(eg)*ed./(real(eg)+ed));
k_x_imag = k0.*imag(eg)./(2*real(eg).^2).*(real(eg)*ed./(real(eg)+ed)).^(3/2);
plot(x,real(e_g),'o',lambda,real(eg))
hold on
plot(x,imag(e_g),'o',lambda,imag(eg))
% 
% 
figure(4)
% h = plot((k_x_real),lambda);
% hold on
% plot((k_x_imag),lambda);
h = plot(real(k_sp),f);
hold on
plot(imag(k_sp),f);
% plot(k0,f)
xlim([0 1e8])
set(gca,'FontSize',10, 'FontName' , 'times new roman','YScale', 'log','XScale', 'log')