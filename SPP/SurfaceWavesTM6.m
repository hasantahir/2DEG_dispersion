function []=SurfaceWavesTM6()
%calculates dispersion relations for electromagnetic waves in an arbitrary homogenous isotropic multilayer film system
%currently written for TM polarization (surface plasmons and some waveguide modes)
%only one mode is solved for, which one depends on starting guess selected by user

clear all
close all

global d
global m
global c
global epsilon
global omega

c=2.998*10^8; %speed of light (m/s)
lambda=1200:-1:500; %wavelength values (nm)
d=[NaN,20,50,20,NaN]; %thickness values in nm (start and end with NaN)
m=length(d); %number of materials
%define refractive index for the different materials and wavelengths (matrix):
n=ones(length(lambda),m);
n(:,2)=sqrt(JCAu(lambda));
n(:,3)=2.24;
n(:,4)=sqrt(JCAu(lambda));

d=d*10^-9; %convert to m
disp('generating imaginary k plot for first wavelength')
k0(1)=2*pi/(lambda(1)*10^-9); %vacuum wavevector
omega=k0(1)*c; %angular frequency
epsilon=n(1,:).^2; %from refractive index to permittivity
stp=400;
detmatrix=zeros(stp,stp);
kxrealmin=0;
kxrealmax=2*10^7;
kxrealstp=(kxrealmax-kxrealmin)/stp;
kximgmin=0;
kximgmax=10^7;
kximgstp=(kximgmax-kximgmin)/stp;
for rlind=1:1:stp
for imind=1:1:stp
kxinput(1)=(kxrealmin+kxrealstp*rlind);
kxinput(2)=(kximgmin+kximgstp*imind);
detmatrix(imind,rlind)=log(MaxwellTM(kxinput));
yrange(imind)=kximgmin+kximgstp*imind;
end
xrange(rlind)=kxrealmin+kxrealstp*rlind;
end
contourf(xrange,yrange,detmatrix,50)
xlabel('real(k) (m-1)')
ylabel('imag(k) (m-1)')
title('click on starting guess')
[x,y]=ginput(1);
kxin=x+1i*y;
disp(['starting guess k = ' num2str(x) ' + ' num2str(y) 'i'])

close all

[k0,kphotons,kx,kz]=dispersion(lambda,n,kxin); %solve dispersion
middle=round(length(lambda)/2); %middle of frequency range
[H,Ex,Ez,z,x]=fields(kz(middle,:),kx(middle),n(middle,:).^2); %calculate field

%plotting:
subplot(1,2,1)
plot(kphotons,k0*c,'b')
hold on
plot(real(kx),k0*c,'r')
axis tight
xlabel('wavevector (m-1)')
ylabel('angular frequency (Hz)')
title('dispersion relation')
subplot(1,2,2)
hold on
plot(lambda,2*pi./real(kphotons)*10^9,'b')
plot(lambda,2*pi./real(kx)*10^9,'r')
axis tight
xlabel('vacuum wavelength (nm)')
ylabel('mode wavelength (nm)')
title('dispersion relation')
figure
plot(lambda,1./imag(kx)*10^6)
axis tight
xlabel('vacuum wavelength (nm)')
ylabel('propagation length (microm)')
title('dissipation')
figure
plot(real(H),z*10^9,'b')
hold on
plot([min(real(H)) max(real(H))],[0 0],'r--')
for a=2:m-1
plot([min(real(H)) max(real(H))],[sum(d(2:a)) sum(d(2:a))]*10^9,'r--') 
end
axis tight
xlabel('field strength (a.u.)')
ylabel('position (nm)')
title('magnetic field');
figure
subplot(2,1,1)
quiver(x*10^9,z*10^9,Ex,Ez,1)
hold on
plot([x(1) x(length(x))],[0 0],'r--')
for a=2:m-1
plot([x(1) x(length(x))],[sum(d(2:a)) sum(d(2:a))],'r--')
end
axis tight
xlabel('position (nm)')
ylabel('position (nm)')
title('E field (vector plot)');
subplot(2,1,2)
surf(x*10^9,z*10^9,real(sqrt(Ex.^2+Ez.^2)))
axis tight
xlabel('position (nm)')
ylabel('position (nm)')
title('E field (intensity plot)');

end

function [k0,kphotons,kx,kz]=dispersion(lambda,n,kxin)

global m
global omega
global epsilon
global c

opt=optimset('TolX',eps,'TolFun',eps,'MaxFunEvals',Inf,'MaxIter',Inf); %settings for fminsearch

%initialize solutions:
kx=zeros(length(lambda),1);
k0=zeros(length(lambda),1);
kphotons=zeros(length(lambda),m);

for a=1:length(lambda) %start solving for each wavelength

k0(a)=2*pi/(lambda(a)*10^-9); %vacuum wavevector
omega=k0(a)*c; %angular frequency
epsilon=n(a,:).^2; %from refractive index to permittivity

%calculate photon lines:
for b=1:m
if real(n(a,b)^2)>0 %check that layer is not metallic
kphotons(a,b)=k0(a)*n(a,b);
else
kphotons(a,b)=NaN; %no wavevector calculated if metallic layer
end
end

if a==1
kxguess=kxin; %first starting guess
else
kxguess=kx(a-1); %starting guess equal to solution for previous wavelength
end

[kxsolved,out,flag]=fminsearch(@MaxwellTM,[real(kxguess),imag(kxguess)],opt);

if flag==1 && a==1
bar=waitbar(0,'dispersion found, solving...');
end
kx(a)=kxsolved(1)+kxsolved(2)*1i;
waitbar(a/length(lambda))

for b=1:m
kz(a,b)=sqrt(epsilon(b)*(omega/c)^2-kx(a)^2);
end

end

close(bar) %get your last drink!

end

function out=MaxwellTM(kxin)

%builds up the matrix describing the field boundary conditions and returns the determinant

global d
global epsilon
global c
global omega
global m

kx=kxin(1)+kxin(2)*1i;

%calculate perpendicular components of wavevectors:
for a=1:m
kz(a)=sqrt(epsilon(a)*(omega/c)^2-kx^2);
end

M=zeros(2*m-2,2*m-2); %build up matrix decsribing boundary conditions
%equations from magnetic field:
M(1,1)=-1;
M(1,2)=exp(-1i*kz(2)*d(2));
M(1,3)=exp(1i*kz(2)*d(2));
for a=2:m-2
M(a,2*a-2)=-1;
M(a,2*a-1)=-1;
M(a,2*a)=exp(-1i*kz(a+1)*d(a+1));
M(a,2*a+1)=exp(1i*kz(a+1)*d(a+1));
end
M(m-1,2*m-4)=-1;
M(m-1,2*m-3)=-1;
M(m-1,2*m-2)=1;
%equations from electric field:
M(m,1)=kz(1)/epsilon(1);
M(m,2)=-kz(2)/epsilon(2)*exp(-1i*kz(2)*d(2));
M(m,3)=kz(2)/epsilon(2)*exp(1i*kz(2)*d(2));
for a=m+1:2*m-3
M(a,2*(a-m))=kz(a-m+1)/epsilon(a-m+1);
M(a,2*(a-m)+1)=-kz(a-m+1)/epsilon(a-m+1);
M(a,2*(a-m)+2)=-kz(a-m+2)/epsilon(a-m+2)*exp(-1i*kz(a-m+2)*d(a-m+2));
M(a,2*(a-m)+3)=kz(a-m+2)/epsilon(a-m+2)*exp(1i*kz(a-m+2)*d(a-m+2));
end
M(2*m-2,2*m-4)=kz(m-1)/epsilon(m-1);
M(2*m-2,2*m-3)=-kz(m-1)/epsilon(m-1);
M(2*m-2,2*m-2)=kz(m)/epsilon(m);

out=det(M);
for a=1:m
out=out/(omega/c); %reduce absolute number in size to improve convergence
end
out=abs(out);

end

function [H,Ex,Ez,z,x]=fields(kz,kx,epsilon)

%solving near field for one wavelength so here kz has m values and kx just one

global d
global m

A=zeros(2*m-2,1);
A(1)=1; %define first value since the overall amplitude is arbitrary

M=zeros(2*m-2,2*m-3); %build up new matrix decsribing boundary conditions
%equations from magnetic field:
M(1,1)=exp(-1i*kz(2)*d(2));
M(1,2)=exp(1i*kz(2)*d(2));
for a=2:m-2
M(a,2*a-3)=-1;
M(a,2*a-2)=-1;
M(a,2*a-1)=exp(-1i*kz(a+1)*d(a+1));
M(a,2*a)=exp(1i*kz(a+1)*d(a+1));
end
M(m-1,2*m-5)=-1;
M(m-1,2*m-4)=-1;
M(m-1,2*m-3)=1;
%equations from electric field:
M(m,1)=-kz(2)/epsilon(2)*exp(-1i*kz(2)*d(2));
M(m,2)=kz(2)/epsilon(2)*exp(1i*kz(2)*d(2));
for a=m+1:2*m-3
M(a,2*(a-m)-1)=kz(a-m+1)/epsilon(a-m+1);
M(a,2*(a-m))=-kz(a-m+1)/epsilon(a-m+1);
M(a,2*(a-m)+1)=-kz(a-m+2)/epsilon(a-m+2)*exp(-1i*kz(a-m+2)*d(a-m+2));
M(a,2*(a-m)+2)=kz(a-m+2)/epsilon(a-m+2)*exp(1i*kz(a-m+2)*d(a-m+2));
end
M(2*m-2,2*m-5)=kz(m-1)/epsilon(m-1);
M(2*m-2,2*m-4)=-kz(m-1)/epsilon(m-1);
M(2*m-2,2*m-3)=kz(m)/epsilon(m);

u=zeros(2*m-2,1);
u(1)=A(1);
u(m)=-kz(1)/epsilon(1)*A(1);

A(2:2*m-2)=(M'*M)\(M'*u); %least squares solution

x=linspace(0,2*pi/real(kx),200); %x range and resolution
z=linspace(-500*10^-9,sum(d(2:m-1))+500*10^-9,200); %z range and resolution

Z=zeros(m);
Z(1)=0;
for a=2:m-1
Z(a)=Z(a-1)+d(a); %the z values of the interfaces
end
Z(m)=Z(m-1);

for a=1:length(z)

layerfound=0;
if z(a)<0 %z=0 defined at the first interface
layerfound=1;
layer=1;
elseif z(a)>sum(d(2:m-1)) %check if z coordinate is at the other side of the layers
layerfound=1;
layer=m;
else
layer=1;
while layerfound==0
layer=layer+1;
if z(a)<sum(d(2:layer))
layerfound=1;
end
end
end

zref=(z(a)-Z(layer)); %z values with reference at each interface

%calcuate magnetic field:
if layer==1
H(a)=A(1)*exp(1i*kz(1)*zref);
elseif layer==m
H(a)=A(2*m-2)*exp(-1i*kz(m)*zref);
else
H(a)=A(2*layer-2)*exp(1i*kz(layer)*zref)+A(2*layer-1)*exp(-1i*kz(layer)*zref);
end

%calculate electric field components:
for b=1:length(x)
if layer==1
Ez(a,b)=-kx/epsilon(1)*A(1)*exp(1i*(kx*x(b)+kz(1)*zref));
Ex(a,b)=kz(1)/epsilon(1)*A(1)*exp(1i*(kx*x(b)+kz(1)*zref));
elseif layer==m
Ez(a,b)=kx/epsilon(m)*A(2*m-2)*exp(-1i*(kx*x(b)+kz(m)*zref));
Ex(a,b)=-kz(m)/epsilon(m)*A(2*m-2)*exp(-1i*(kx*x(b)+kz(m)*zref));
else
Ez(a,b)=-kx/epsilon(layer)*A(2*layer-2)*exp(1i*(kx*x(b)+kz(layer)*zref))+kx/epsilon(layer)*A(2*layer-1)*exp(-1i*(kx*x(b)+kz(layer)*zref));    
Ex(a,b)=kz(layer)/epsilon(layer)*A(2*layer-2)*exp(1i*(kx*x(b)+kz(layer)*zref))-kz(layer)/epsilon(layer)*A(2*layer-1)*exp(-1i*(kx*x(b)+kz(layer)*zref));
end
end

end

end

function epsilon=JCAu(lambda)

%analytical formula for gold with wavelength in nm, fits J&C data:
epsiloninf=1.54;
lambdap=143;
gammap=14500;
A1=1.27;
lambda1=470;
phi1=-pi/4;
gamma1=1900;
A2=1.1;
lambda2=325;
phi2=-pi/4;
gamma2=1060;

for a=1:length(lambda)
epsilon(a)=epsiloninf-1/(lambdap^2*(1/lambda(a)^2+1i/(gammap*lambda(a))))...
+A1/lambda1*(exp(phi1*1i)/(1/lambda1-1/lambda(a)-1i/gamma1)+exp(-phi1*1i)/(1/lambda1+1/lambda(a)+1i/gamma1))...
+A2/lambda2*(exp(phi2*1i)/(1/lambda2-1/lambda(a)-1i/gamma2)+exp(-phi2*1i)/(1/lambda2+1/lambda(a)+1i/gamma2));
end

end