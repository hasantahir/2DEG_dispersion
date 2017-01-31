% %% page # 693 Evaluation of Definite Integrals by Use of the 
% %   Residue Theorem
% % Boas - Mathematical Methods in the Physics
% %
% %
% %
%  integral (r^(p-1)/(r+1)),0,inf)
clear all; close all
R = 2; 
a = .5;
p = .2; % 0 < p < 1
sumc = 0;
sum1 = 0;
nstep = 6000;
ans_p = -2*pi*1i*exp(1i*pi*p);
alpha = asin(a/R);
thetai = 0 + alpha;
thetaf = 2*pi - alpha;
dtheta = (thetaf - thetai)/nstep;

% Integration around the circle
for i = 1 : nstep
    theta = thetai + dtheta*i;
    z = R*exp(1i*theta);
    dz = R*exp(1i*theta)*1i*dtheta;
    f = @(z) abs(z).^((p-1).*exp(1i*theta*(p-1)))...
    ./(1+abs(z).^(exp(1i*theta)));
    sumc = sumc+(dz/2.)*( feval(f,z) + feval(f,z-dz) );
end
% Vertical Integration
dz = 2*1i*sqrt(R^2 - a^2)/nstep;
for i = 1 : nstep
    z = -a -1i*sqrt(R^2-a^2)+ dz*i;
    if imag(z) > 0
        theta = 0 + atan(a/abs(imag(z)));
    elseif imag(z) < 0
        theta = 2*pi - atan(a/abs(imag(z)));
    end
     sum1 = sum1+(dz/2.)*( feval(f,z) + feval(f,z-dz) );
end
% disp(z)
partR = real(sumc)+ real(sum1)
partI = imag(sumc)+ imag(sum1)
syms x
X = @(x) x.^(p-1)./(x+1);
integral(X,0,inf)

% % c0 = chebfun('1.5*exp(1i*pi*s)',[0.51 2.49]);          % big circle
% c0 = chebfun(@(s) 1.5*exp(1i*pi*s), [ .51 2.49]);
% % c1 = chebfun('1+.2*exp(-1i*pi*s)',[-0.93 0.93]);       % right circle
% c1 = chebfun(@(s) (1+.2*exp(-1i*pi*s)),[-0.93 0.93]);
% c2 = -c1;                                              % left circle
% p1 = c0(0.51); p2 = c0(2.49);
% p3 = real(c0(2.49)) + 1i*imag(c1(-0.93));
% p4 = c1(-0.93); p5 = c1(0.93);                         % corner points
% p6 = c2(-0.93); p7 = c2(0.93);
% p8 = real(c0(0.51)) + 1i*imag(c2(0.93));
% s = chebfun('s',[0 1]);
% z = join( c0, p2+s*(p3-p2), p3+s*(p4-p3), c1, ...      % the contour
%        p5+s*(p6-p5), c2, p7+s*(p8-p7), p8+s*(p1-p8) );
% plot(z,'k','linewidth',1.6), ylim([-1.8 1.8])
% hold on, plot([-1 1],[0 0],'.r','markersize',10), hold off
% axis equal, title('Ablowitz-Fokas contour','fontsize',14)
% ff = @(z) (.5i/pi)*(z.^2-1).^(1/2).*(-1).^(real(z)>0)./(1+z.^2);
% % I = sum(ff(z).*diff(z))