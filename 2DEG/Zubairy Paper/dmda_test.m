clear all; close all
lambda = linspace(12e2, 60e5, 1e3); % wavelength in lambda

c = 3e8;
% k0 = 2*pi/la0; 
a = 40; % la0,a in units of nm
es = 12.8; 
ec = 10.7; 
ef = -.113 - .001*1i;%(0.0657-4*1i)^2;
tol = 1e-10; % error tolerance

be0 = zeros(size(lambda));
E0 = zeros(size(lambda));
N0 = zeros(size(lambda));
for i = 1 : length(lambda)
    k0 = 2*pi/lambda(i); 
    [be0(i),E0(i),N0(i)] = dmda(lambda(i),ef,ec,es,a,0,tol); % LRSP
end
plot(real(be0), c./(lambda*1e-9),'o');
hold on
plot(2*pic./(lambda*1e-9), c./(lambda*1e-9))
% plot(imag(be0), c./lambda);
set(gca,'Xscale','Log','Yscale','Log')
