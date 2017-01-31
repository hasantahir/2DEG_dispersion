clear;clf
q = 1;
format long
while q>0
    xo = input('x coord center of box=');
    yo = input('y coord center of box=');
    dx = input('half width of box, x direction=');
    dy = input('half width of box, y direction=');
    %below use odd numbers for nx and ny to include center of box
    %in your calculations
    nx=input('number of x divisions');
    ny=nx; %this is numberof y divisions
    % the two lines below give the real and imaginary parts of c.
    tic
    cr=linspace(xo-dx,xo+dx,nx);
    ci=linspace(yo-dy,yo+dy,ny);
    [Cr,Ci]=meshgrid(cr,ci);
    c=Cr+1i*Ci;%this creates a grid of complex numbers for c
    nmax=1000;% we use 1000 iterations of the recursion relation
    j=1;
    z=zeros(size(c));%this starts off z, for each value of c, at the value zero
    while j<=nmax;%iterates the expression below nmax times;
        z=z.*z+c;
        j=j+1;
    end
    ck=abs(z)<=2;%puts a symbol of 1 in the matrix where |z|<=2 ; otherwise puts %0
    dk=1.0*ck;%converts symbolic elements to numerical in above matrix.
    [rows,cols,vals] = find(dk);%this finds the nonzero elements in dk
    for k=1:length(rows)
        locations=c(rows(k),cols(k));
        figure(1)
        plot(locations,'k.'); hold on
        %the above is optional and will plot the points
        %you just found as dots
        %the "locations" are points in the complex plane lying in the
        % mandelbrot set
    end
    toc
    q = input('negative q to stop');
end