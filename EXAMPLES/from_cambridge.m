% p21: square root branch cut; 04/02/2008 SGLS
clf;
x = -4:.01:4; y = x;
[x,y] = meshgrid(x,y); z = x + 1i*y;
beta = 1 + .1i;
f = sqrt(z.^2-beta.^2);
% f = log(z);
% f = z.^2;
surfc(x,y,imag(f)); shading interp
xlabel('Re z'); ylabel('Im z'); zlabel('Im $\sqrt {z}$','Interpreter','latex')
colormap viridis
view([0 90])
% contour(imag(f));
