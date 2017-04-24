clear all;close all;clc

% Read an image
X = imread('new.png');
z = rgb2gray(X);
z = double(z);

% define coordinaties of z
x = (1:length(z))*2e-9;
y = (1:length(z))*2e-9;

% Create a numerical space
len = 5e1;
f = 25e12; % free-space frequency
c = 3e8;
lambda0 = c/f;
k0 = 2*pi/lambda0;

k = linspace(-k0, k0, len);
kx = k;
ky = k;


fZ = zeros(length(kx), length(ky));
ifz = zeros(size(z));


% Compute Fourier Transform
%% Display Progress Bar
h = waitbar(0,'Calculating DFT...');
steps = length (k);
for u = 1 : length(k)
    for v = 1 : length(k)
        for i = 1 : length(z)
            for j = 1 : length(z)
                
                dft_exp = -2 * 1i * pi*(kx(u) * x(i) + ky(v) * y(j));
                fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	% D.F.T
            end
        end
    end
    waitbar(u / steps);
end
set(h, 'HandleVisibility', 'on')
% I.D.F.T.
h = waitbar(0,'Calculating inverse DFT...');
steps = length (z);
for i = 1 : length(z)
    for j = 1 : length(z)
        for u = 1 : length(k)
            for v = 1 : length(k)
    
                idft_exp = +2 * 1i * pi*(kx(u) * x(i) + ky(v) * y(j));
                ifz(i,j) = ifz(i,j) + fZ(u,v) * exp(idft_exp);	%%%!incident field dft
            end
        end
    end
    waitbar(i / steps);
end
set(h, 'HandleVisibility', 'on')
close all


figure(1)
surf(z)
shading interp; colormap(jet(3));
view([0 90])
axis tight

figure(2)
sfZ = fftshift(fZ);
surf(abs(sfZ))
shading interp; colormap(jet);
view([0 90])
axis tight

figure(3)
surf(abs(ifz)) ; 
shading interp; colormap(jet);
view([0 90])
axis tight
        