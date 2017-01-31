% % p21: square root branch cut; 04/02/2008 SGLS
clf; close all;clear all
x = -4:.001:4; y = x;
[x,y] = meshgrid(x,y); z = x + 1i*y;
beta = 1 + .5i;
% f = 1./sqrt(1-z.^2)./(1+z.^2);
% f = sqrt(z.^2-sqrt(z.^2-beta.^2));
f = sqrt(z.^2+beta.^2) + sqrt(z.^2+real(beta).^2);
% f = log(z);
% f = z.^2;
% f = z.^2;
figure(1)
surf(x,y,real(f)); shading interp
xlabel('Re z'); ylabel('Im z'); zlabel('Im $\sqrt {z}$','Interpreter','latex')
colormap viridis
view([0 90])
figure(2)
% contour(x,y,imag(f),[1],...
%     'LineWidth' , 1.4,...
%     'Color','black');
% contour(x,y,imag(f),[1],...
%     'LineWidth' , 1.4,...
%     'Color','black');
% p =.4;
% f = @(z) 1./sqrt(1-z.^2)./sqrt(1+z.^2);
% % f = @(z) z.^p./(1 + sqrt(2)*z + z.^2);
% X = quadgk(f,0,inf,...
%     'MaxIntervalCount',20e5,...
%     'RelTol',0,...
%     'AbsTol',1e-10)
% Y = sqrt(2)*pi*sin(p*pi/4)/sin(p*pi);
% X-Y
% YY = 1.3110 - 1.3110i;
% p32: iterations and fractals; 04/07/2008 SGLS
% x = -1.2:0.0005:1.2; y = x;
% [x,y] = meshgrid(x,y); z = x + i*y;
% c = 0.285 + 0.01*i;
% for j = 1:length(x);
%     for k = 1:length(y);
%         zn = z(j,k); n = 0;
%         while (abs(zn)<10 & n<100)
%             zn = zn.^2 + c; n = n + 1;
%         end
%         f(j,k) = n;
%     end
% end
% pcolor(x,y,f); shading flat; axis square; colorbar
% xlabel('Re z'); ylabel('Im z');
