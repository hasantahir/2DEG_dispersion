close all;clear all
%% First example from the Numerical resulsts section of Michalski's review paper
h = 5 ;
eps_r = 10 - 1i*180;
c = 2.99792458e8;
xmu = 4*pi*1e-7;
eps_0 = 8.854187817e-12;
lambda_0 = 300;
f_0 = 1e6;
omega = 2*pi*f_0;
sigma = 10e-3;

k_1 = 2*pi/lambda_0;
k_2 = omega*sqrt(eps_r*eps_0*xmu);
rho = 50;
% R_p = -1i*(eps_r*sqrt(eps_r))/(eps_r^2 - 1)*exp(1i*k_1*h/(sqrt(eps_r + 1)))
% Pole location
% x = 1 - (real(eps_r) + 3/4) / (2*imag(eps_r)^2) - 1i/(2*imag(eps_r))
real_kp_by_k1 = .4: .01 : 1.2;
imag_kp_by_k1 = -.4: .01 : .4;
[x,y] = meshgrid(real_kp_by_k1, imag_kp_by_k1);
k_p = k_1*(x + 1i*y);
k_z1 = sqrt(k_1^2 - k_p.^2);
k_z2 = sqrt(k_2^2 - k_p.^2);
f = -1i*eps_r./(eps_r*k_z1 + k_z2).*exp(-1i*k_z1*h);
J = besselj(0,rho*k_p);
S_0 = f.*J.*k_p;
axes1 = axes('Parent',gcf);
hold(axes1,'on');
surfc(x,y,log(abs(S_0)))
% Create axes
hold on
shading interp
colormap viridis
view(axes1,[54.1 36.4]);
grid(axes1,'on');
axis(axes1,'tight');
% Set the remaining axes properties
set(axes1,'CameraViewAngle',9.12883787926109,'DataAspectRatio',...
    [1 1 4.17798655114009],'Projection','perspective');
alpha .5
plot3(real_kp_by_k1,y(41,:),log(abs(S_0(41,:))),'LineWidth',1.2, 'Color','black')
contour(x,y,log(abs(S_0)))

% kz1 = @(kp) sqrt(k_1^2 - kp^2);
% kz2 = @(kp) sqrt(k_2^2 - kp^2);
% f = @(kp) 1i*eps_r./(eps_r*k_z1 + k_z2).*exp(-1i*k_z1*h);
% J = @(kp) besselj(0,rho*kp)
% syms r
% rho = linspace(0,2,21);
% hankel3 = @(r,rho) exp(-r).*besselj(0,2*pi*rho*r);
% 
% g = zeros(size(rho));
% for k=1:length(rho)
%     g(k) = 2*pi*quad(hankel3,0,5,[],[],rho(k));
% end
% 
% u = linspace(0,2,101);
% v = 2*pi./sqrt(4*(pi*u).^2+1);
% plot(u,v,rho,g,'ko');
% xlabel('radial freq \rho');
% ylabel('value');
% g = zeros(size(rho));
% for k=1:length(rho)
% g(k) = quad(hankel3,0,5,rho(k));
% end
