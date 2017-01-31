% c0 = chebfun('1.5*exp(1i*pi*s)',[0.51 2.49]);          % big circle
% % c0 is defined by the expression on 0.51 to 2.49
% c1 = chebfun('1+.2*exp(-1i*pi*s)',[-0.93 0.93]);       % right circle
% c2 = -c1;                                              % left circle
% p1 = c0(0.51);
% p2 = c0(2.49);
% p3 = real(c0(2.49)) + 1i*imag(c1(-0.93));
% p4 = c1(-0.93);
% p5 = c1(0.93);                        % corner points
% p6 = c2(-0.93);
% p7 = c2(0.93);
% p8 = real(c0(0.51)) + 1i*imag(c2(0.93));
% s = chebfun('s',[0 1]);
% z = join( c0, p2+s*(p3-p2), p3+s*(p4-p3), c1, ...      % the contour
%        p5+s*(p6-p5), c2, p7+s*(p8-p7), p8+s*(p1-p8) );
% plot(z,'k','linewidth',1.6), ylim([-1.8 1.8])
% hold on, plot([-1 1],[0 0],'.r','markersize',10), hold off
% axis equal, title('Ablowitz-Fokas contour','fontsize',14)
% 
% ff = @(z) (.5i/pi)*(z.^2-1).^(1/2).*(-1).^(real(z)>0)./(1+z.^2);
% figure(3)
% c = (-1).^(real(z)>0)
% plot(c)
% I = sum(ff(z).*diff(z))
clf; clear
% %% Example 4.3.4 Ablowitz and Fokas
a = 1;
k = .5;
CR = chebfun('1.5*exp(1i*s)',[0.01*pi 2*pi-.01]);          % big circle
% cR is defined by the expression on 0.51 to 2.49
C1 = chebfun('.2*exp(-1i*s)',[.1*pi 2*pi-.1]);         % small circle
P1 = CR(0.01*pi); % right-most point on the big circle in first quadrant
P2 = CR(2*pi-.01); % right-most point on the big circle in fourth quadrant
P3 = C1(.1*pi); % right-most point on the small circle in first quadrant
P4 = C1(2*pi-.1); % right-most point on the small circle in fourth quadrant
% 
% % p4 = c1(-0.93); p5 = c1(0.93);                         % corner points
% % p6 = c2(-0.93); p7 = c2(0.93);
% % p8 = real(c0(0.51)) + 1i*imag(c2(0.93));
% figure(2)
s = chebfun('s',[0 1]);
z = join(CR,P2-s*(P2-P3), C1, P4+s*(P1-P4));
plot(z,'k','linewidth',1.6), ylim([-1.8 1.8])
hold on, plot(0,1,'xr','markersize',10), hold off;axis equal
hold on, plot(0,-1,'xr','markersize',10), hold off;axis equal
hold on, plot(0,0,'ob','markersize',4), hold off;axis equal
ff1 = @(z) z.^(-k)./((a+z.^2).^2).*(-1).^(real(z)>0);
% % % 'Waypoints'
% % % where the factor involving real(z)
% % % appears in order to avoid inappropriate jumps
% % % of branch when zz crosses the negative imaginary axis.
 I1 = sum(ff1(z).*diff(z))
%%
% x = 10e10:1e17:1e20;
% syms t
% I = zeros(size(x));
% for i = 1 : length(x)
%     J = @(t) 1./pi*cos(x(i).*sin(t));
%     I(i) = integral(J,0,pi);
% end
    
    
    

