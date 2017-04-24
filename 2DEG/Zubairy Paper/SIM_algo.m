% Create a numerical space
len = 1e3;
x = linspace (0, 1e-6, len);
y = linspace (0,1e-6, len);

[z] = meshgrid(x,y);
fZ = zeros(size(z));
% Create objects as small beads randomly

% 1st object
for i = 1 : len
    for j = 1 : len
        if x(i) < 500e-9 && x(i) >  450e-9 ...
                && y(j) < 500e-9 && y(j) >  450e-9
            z(i,j) = 1;
        
        else
            z(i,j) = 0;
        end
    end
end

x = linspace (0, .5e-6, len);
y = linspace (0,.5e-6, len);

[z] = meshgrid(x,y);
fZ = zeros(size(z));
% Create objects as small beads randomly

% 1st object
for i = 1 : len
    for j = 1 : len
        if x(i) < 260e-9 && x(i) >  240e-9 ...
                && y(j) < 260e-9 && y(j) >  240e-9
            z(i,j) = 1;
        
        else
            z(i,j) = 0;
        end
    end
end



for u = 1 : len
    for v = 1 : len
        for x = 1 : len
            for y = 1 : len
    
                dft_exp = -2 * 1i * pi * (u * x / len + v * y / len);
                fZ(u,v) = fZ(u,v) + z(x,y) * exp(dft_exp);	%%%!incident field dft
            end
        end
    end
end

I = mat2gray(z, [0 1]);
imshow(I)
cmap = jet(8);
colormap(cmap);


            