%  clear all; a = 1000;
% t = 2;
%  % C PATH
% f = @(z) besselh(0,1,z).*exp(-1j*z*t);
% z = @(R) R - sqrt(a^2 - R.^2)*1j;
% dzdr = @(R) 1 - 1/2*(a^2 - R.^2).^(-1/2).*(-2.*R)*1j;
% g_C = integral(@(R) f(z(R)).*dzdr(R), a, -a)
% g_C_2 = integral(@(R) real(f(z(R))).*dzdr(R), a, -a) + 1j*integral(@(R) imag(f(z(R))).*dzdr(R), a, -a)
% g_C_3 = integral(@(R) real(f(z(R)).*dzdr(R)), a, -a) + 1j*integral(@(R) imag(f(z(R)).*dzdr(R)), a, -a) 
% % I PATH
% f = @(z) besselh(0,1,z).*exp(-1j*z*t);
% z = @(R) R;
% dzdr = @(R) 1;
% g_I = integral(@(R) f(z(R)).*dzdr(R), a, -a)
% % I PATH
% f = @(z) besselh(0,1,z).*exp(-1j*z*t);
% z = @(R) R - 1e-20j; % under branch cut
% dzdr = @(R) 1;
% g_I = integral(@(R) f(z(R)).*dzdr(R), a, -a)
% % C PATH
% f = @(z) besselh(0,1,z,1).*exp(-1j*z*(t-1));
% z = @(R) R - sqrt(a^2 - R.^2)*1j;
% dzdr = @(R) 1 - 1/2*(a^2 - R.^2).^(-1/2).*(-2.*R)*1j;
% g_C = integral(@(R) f(z(R)).*dzdr(R), a, -a)

% %% Integrate function from 0 to infty
% % f(x) = x^k         -1< k < 3
% %       ______     
% %        (x^2+1)^2
% % syms x
% k = 1;
% z = @(R,theta) R.*exp(1i*theta); 
% f = @(R,theta) (R.*exp(1i*theta)).^k./((R.*exp(1i*theta)).^2 + 1).^2;
% I1 = integral(@(R) f(R,eps),0+eps,inf,'RelTol',1e-12,'AbsTol',1e-12);
% I2 = integral(@(R) f(R,2*pi - eps),inf,0+eps,'RelTol',1e-12,'AbsTol',1e-12);
% I = I1 + I2 
% pi/4*(1-k)/(cos(pi*k/2))






clear all



% c0 = chebfun('1.5*exp(1i*pi*s)',[0.51 2.49]);          % big circle
c0 = chebfun(@(s) 1.5*exp(1i*pi*s), [ eps 2*pi-eps]);
% c1 = chebfun('1+.2*exp(-1i*pi*s)',[-0.93 0.93]);       % right circle
c1 = chebfun(@(s) (0.2*exp(1i*pi*s)),[2*pi-eps eps]);
p1 = c0(eps); p2 = c0(2*pi-eps);
p3 = real(c0(eps)) + 1i*imag(c1(-0.93));
p4 = c1(-0.93); p5 = c1(0.93);                         % corner points
p6 = c2(-0.93); p7 = c2(0.93);
p8 = real(c0(0.51)) + 1i*imag(c2(0.93));
s = chebfun('s',[0 1]);
z = join( c0, p2+s*(p3-p2), p3+s*(p4-p3), c1, ...      % the contour
       p5+s*(p6-p5), c2, p7+s*(p8-p7), p8+s*(p1-p8) );
plot(z,'k','linewidth',1.6), ylim([-1.8 1.8])
hold on, plot([-1 1],[0 0],'.r','markersize',10), hold off
axis equal, title('Ablowitz-Fokas contour','fontsize',14)
ff = @(z) (.5i/pi)*(z.^2-1).^(1/2).*(-1).^(real(z)>0)./(1+z.^2);
I = sum(ff(z).*diff(z))
