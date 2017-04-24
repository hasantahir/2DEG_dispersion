clear all;close all;clc

% Read an image
X = imread('new.png');
z = rgb2gray(X);
z = double(z);
z = padarray(z, [2,2],'replicate','post');

% define coordinaties of z
x = (1:length(z))*1e-9;
y = (1:length(z))*1e-9;

% Create a numerical space
len = 3e1;
f = 25e12; % free-space frequency
c = 3e8;


% Frequency domain
lenk = 1e2;
lambda0 = c/f; % free-space wavelength
k0 = 2*pi/lambda0;

lambdap = .15e-6; % plasmonic wavelength
kp = 2*pi/lambdap;

k = linspace(0, k0, lenk);
kx = k;
ky = k;

kp = linspace(0, kp, lenk);
kxp = kp;
kyp = kp;

fZ = zeros(length(kx), length(ky));
ifz = zeros(size(z));

fZp = zeros(length(kxp), length(kyp));
ifzp = zeros(size(z));

% Dimensions
M = len * max(x);
N = len * max(y);

% D.F.T.
for u = 1 : lenk
    for v = 1 : lenk
        for i = 1 : len
            for j = 1 : len
    
                dft_exp = -2 * 1i * pi * (kx(u) * x(i) / M + ky(v) * y(j) / N);
                fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	%%%!incident field dft
            end
        end
    end
end


% I.D.F.T.
for i = 1 : len
    for j = 1 : len
        for u = 1 : lenk
            for v = 1 : lenk
    
                idft_exp = 2 * 1i * pi * (u * i / M + v * j / N);
                ifz(i,j) = ifz(i,j) + fZ(u,v) * exp(idft_exp);	%%%!incident field dft
            end
        end
    end
end

figure(1)
surf(z)
shading interp; colormap(jet(3));
view([0 90])
axis tight

figure(2)
sfZ = fftshift(fZ);
surf(abs(fZ))
shading interp; colormap(jet);
view([0 90])
axis tight

figure(3)
surf(1/(M*N)*abs(ifz)) ; 
shading interp; colormap(jet);
view([0 90])
axis tight