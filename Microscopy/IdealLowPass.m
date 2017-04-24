% Program to make Moire patterns.
% First Moire pattern is from interference (multipication) to two linear sine wave ripples.
% This is equivalent to the Fourier diffraction pattern of two infinitely long slits.
% The second Moire pattern is from interference (multipication) to two Sombrero function ripples.
% This is equivalent to the Fourier diffraction pattern of two infinitely small circular apertures.

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
% workspace;  % Make sure the workspace panel is showing.
% format long g;
% format compact;
fontSize = 12;
rows = 256;
columns = 256;

% Make an image with rippled lines in it.

rowVector = (1 : rows)';
period = 30; % 20 rows
cosVector = 10*(cos(2 * pi * rowVector / period));
ripplesImage = repmat(cosVector, [1, columns]);
imshow(ripplesImage,[])
axis square
colormap(brewermap([],'Spectral'))


sss = double(imread('fine_lines.png'));
subplot(2, 3, 1);
imshow(ripplesImage,[])
% imshow(sss, []);
colormap(brewermap([],'Spectral'))
axis on;
title('Ripple image', 'FontSize', fontSize);



% Make an image with tilted rippled lines in it.
ripplesImage2 = imrotate(ripplesImage, 10.6, 'crop'); % Rotate first image.
subplot(2, 3, 2);
imshow(ripplesImage2, []);
colormap(brewermap([],'Spectral'))
axis on;
title('Tilted Ripples', 'FontSize', fontSize);

% Multiply the ripple images together to get a Moire pattern.
grayImage = ripplesImage .* (ripplesImage2);
subplot(2, 3, 3);
imshow(grayImage, []);
colormap(brewermap([],'Spectral'))
axis on;
title('Moire Image', 'FontSize', fontSize);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
drawnow;

y1 = fft2(ripplesImage,2048,2048);
y1 = fftshift(y1);
subplot(2, 3, 4);
imshow((abs(y1)), [])
colormap(brewermap([],'Spectral'))
% axis on;
title('FT of Ripple image', 'FontSize', fontSize);

y2 = fft2(ripplesImage2,2048,2048);
y2 = fftshift(y2);
subplot(2, 3, 5);
imshow((abs(y2)), [])
colormap(brewermap([],'Spectral'))
% axis on;
title('FT of Rotated image', 'FontSize', fontSize);

y3 = fft2(grayImage,2048,2048);
y3 = fftshift(y3);
subplot(2, 3, 6);
imshow((abs(y3)), [])
colormap(brewermap([],'Spectral'))
% axis on;
title('FT of Moire image', 'FontSize', fontSize);