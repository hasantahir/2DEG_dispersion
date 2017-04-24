clear all;close all;clc

% Read an image
X = imread('new.png');
figure('Name','RGB Image')
imshow(X)
axis image; axis on
x = rgb2gray(X);
x = double(x);

% Create a numerical space
len = 1e2;
f = 25e12; % free-space frequency
c = 3e8;
k0 = 2*pi*f/c;

k = linspace(0, k0, len);


i = linspace (0, 1e-6, len);
j = linspace (0,1e-6, len);

[z] = meshgrid(i,j);
fZ = zeros(size(z));
ifz = zeros(size(z));


% Create objects as small beads randomly

% 1st object
for i = 1 : len
    for j = 1 : len
        if i(i) < 260e-9 && i(i) >  240e-9 ...
                && j(j) < 260e-9 && j(j) >  240e-9
            z(i,j) = 1;
        
        else
            z(i,j) = 0;
        end
    end
end





% D.F.T.
for u = 1 : len
    for v = 1 : len
        for i = 1 : len
            for j = 1 : len
    
                dft_exp = -2 * 1i * pi * (k(u) * x(i) / len + v * j / len);
                fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	%%%!incident field dft
            end
        end
    end
end

% I.D.F.T.
for i = 1 : len
    for j = 1 : len
        for u = 1 : len
            for v = 1 : len
    
                idft_exp = 2 * 1i * pi * (u * i / len + v * j / len);
                ifz(i,j) = ifz(i,j) + fZ(u,v) * exp(idft_exp);	%%%!incident field dft
            end
        end
    end
end

I = mat2gray(z, [0 1]);
imshow(I)
cmap = jet(8);
colormap(cmap)
axis auto
sfZ = fftshift(fZ);
surf(abs(sfZ));
surf(abs(ifz)) ;           