clear all;close all;clc

X = imread('new.png');
figure('Name','RGB Image')
imshow(X)
axis image; axis on
x = rgb2gray(X);
x = double(x);
y = fft2(x,2000,2000);
y = fftshift(y);


% D.F.T.
% for u = 1 : len
%     for v = 1 : len
%         for i = 1 : len
%             for j = 1 : len
%     
%                 dft_exp = -2 * 1i * pi * (k(u) * x(i) / len + v * j / len);
%                 fZ(u,v) = fZ(u,v) + z(i,j) * exp(dft_exp);	%%%!incident field dft
%             end
%         end
%     end
% end

imshow(log(abs(y)), [])
colormap jet
yy = ifftshift(y);
xx = ifft2(yy);
xx = imshow(abs(xx));
colormap jet(3)

% [f1,f2] = freqspace(1024,'meshgrid');
% Hd = zeros(1024,1024); d = sqrt(f1.^2 + f2.^2) < .2;
% Hd(d) = 1;
% mesh(f1,f2,Hd)
% 
% YY = y.*Hd;
% yyy = ifftshift(YY);
% xx = ifft2(yyy);
% xx = imshow(real(xx));
