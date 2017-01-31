[x, y] = meshgrid(-20:.1:10,-20:.1:10);
z = x + 1i*y;
z = besselh(0,1,z);
figure(1)
surf(x,y,real(z))
shading interp
zlim([-1 1])
view([145 60])
caxis([-100 100])
figure(2)
surf(x,y,imag(z))
shading interp
zlim([-1 1])
view([145 60])
caxis([-100 100])
fun=@(z)(z+1)./(z.^3-2*z.^2);     
g=@(theta)2+cos(theta)+1i*sin(theta);
gprime=@(theta)-sin(theta)+1i*cos(theta);
q1 = integral(@(t) fun(g(t)).*gprime(t),0,2*pi)