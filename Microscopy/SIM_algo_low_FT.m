clear all;close all;clc

% Read an image
X = imread('new.png');
z = rgb2gray(X);
z = double(z);
z = padarray(z, [2,2],'replicate','post');

% define coordinaties of z
x = 1:length(z);
y = 1:length(z);

% Create a numerical space
len = 3e1;
f = 25e12; % free-space frequency
c = 3e8;

lambda0 = c/f; % free-space wavelength
k0 = 2*pi/lambda0;

lambdap = .15e-6; % plasmonic wavelength
kp = 2*pi/lambdap;

k = linspace(-k0, k0, len);
kx = k;
ky = k;

kp = linspace(-kp, kp, len);
kxp = kp;
kyp = kp;

fZ = zeros(length(kx), length(ky));
ifz = zeros(size(z));

fZp = zeros(length(kxp), length(kyp));
ifzp = zeros(size(z));


% % Compute Fourier Transform
% for i = 1 : length(z)
%     for j = 1 : length(z)
%         for u = 1 : length(k)
%             for v = 1 : length(k)  
%                 dft_exp = -2 * 1i * pi*(kx(u) * x(i)*1e-9 + ky(v) * y(j)*1e-9);
%                 fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	% D.F.T
%             end
%         end
%     end
% end

% Compute Fourier Transform
%% Display Progress Bar
h = waitbar(0,'Calculating DFT...');
steps = length (k);
for u = 1 : length(k)
    for v = 1 : length(k)
        for i = 1 : length(z)
            for j = 1 : length(z)
                
                dft_exp = -2 * 1i * pi*(kx(u) * x(i)*2e-9 + ky(v) * y(j)*2e-9);
                fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	% D.F.T
            end
        end
    end
    waitbar(u / steps);
end

% I.D.F.T.
h = waitbar(0,'Calculating inverse DFT...');
steps = length (z);
for i = 1 : length(z)
    for j = 1 : length(z)
        for u = 1 : length(k)
            for v = 1 : length(k)
    
                idft_exp = 2 * 1i * pi*(kx(u) * x(i)*2e-9 + ky(v) * y(j)*2e-9);
                ifz(i,j) = ifz(i,j) + fZ(u,v) * exp(idft_exp);	%%%!incident field dft
            end
        end
    end
    waitbar(i / steps);
end

%% Plasmonic DFT
% h = waitbar(0,'Calculating DFT...');
% steps = length (kp);
% for u = 1 : length(kp)
%     for v = 1 : length(kp)
%         for i = 1 : length(z)
%             for j = 1 : length(z)
%                 
%                 dft_exp = -2 * 1i * pi*(kxp(u) * x(i)*2e-9 + kyp(v) * y(j)*2e-9);
%                 fZp(u,v) = fZp(u,v) + z(i,j) * exp(dft_exp);	% D.F.T
%             end
%         end
%     end
%     waitbar(u / steps);
% end
% 
% % Plasmonic I.D.F.T.
% h = waitbar(0,'Calculating inverse DFT...');
% steps = length (z);
% for i = 1 : length(z)
%     for j = 1 : length(z)
%         for u = 1 : length(kp)
%             for v = 1 : length(kp)
%     
%                 idft_exp = 2 * 1i * pi*(kx(u) * x(i)*2e-9 + ky(v) * y(j)*2e-9);
%                 ifz(i,j) = ifz(i,j) + fZ(u,v) * exp(idft_exp);	%%%!incident field dft
%             end
%         end
%     end
%     waitbar(i / steps);
% end
close(h);

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
surf(abs(ifz)) ; 
shading interp; colormap(jet);
view([0 90])
axis tight
        