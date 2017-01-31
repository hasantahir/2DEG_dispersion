% close all;clf;
% % xmin=-0.5; xmax=0.5; ymin=-0.5; ymax=0.5;
% % xres = 400; yres = 400;
% % x = linspace(xmin,xmax,xres);
% % y = linspace(ymin,ymax,yres);
% % z = ones(yres,1)*x + 1i*y'*ones(1,xres);
% % f = exp(1./z);
% % p = surf(real(z),imag(z),0*f,angle(-f));
% % set(p,'EdgeColor','none');
% % caxis([-pi,pi]), colormap hsv(600)
% % view([0 90])
% 
n = 2;
m = 50;
r = (0:m)'/m;
theta = pi*(-n*m:n*m)/m;
z = r * exp(1i*theta);
% s = r.^(1/n) * exp(1i*theta/n);
s = log(r) * exp(1i*theta/n);
% s = (z.^2 - 1).^(.5)./(1+z.^2);
surfc(real(z),imag(z),real(s),'EdgeColor','interp','FaceLighting','gouraud');
shading flat % Avoid jittered shading
caxis([-1,1]) % set the colorbar range from -1 to 1
set(gcf,'Color','white'); % Set background color to white
title('Riemann Surface of $f(z) = \sqrt{z}$','Interpreter','latex')
set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
ax = gca;
xlabel('$\Re(z)$','Interpreter','latex');
ylabel('$\Im(z)$','Interpreter','latex');
% material dull; % Set reflectivity of the surface to dull
colormap viridis
% alpha .7
% view([ -120 30])
% zlim([ 0 10])
% axis equal
% matlab2tikz('filename',sprintf('Square_root_Riemann.tex'))
% 
% [x,w]=lgwt(5,0,pi);
% F = 0;
% for i = 1 : length(x)
%     F = F + sum(sin(x(i))*w(i));
% end
% F
% syms x;f = @(x) sin(x);
% f1 = quadgk(f,0,pi)