%  zmap.m
%
%  Contour magnitude and phase of an analytic function f(z) on
%  the complex plane.  f(z) is given by a matlab M-file, fun.m 
%  which the user must provide.  An example file is provided, which
%  can be edited to the user's specifications.
%
%  There are several viewing options: see comments in
%    solid.m  zview.m  turn.m  zooom.m  pullback.m  pan.m
%
%  Values of f(z) are saved in the complex array  FZ 
%
%  Integration in the complex plane on a path drawn on the map
%  can be performed by the program  integrate.m
%

%  Center and side of window set here
if (exist('center') == 0)
  side=[6.3, 6.3]
  center=[0 0]
  nside=81;
end


%  Set up a grid in the complex plane
ds = max([side(1), side(2)]) / (nside-1);
nx = 1+side(1)/ds;  ny = 1+side(2)/ds;
xx = center(1) - side(1)/2 + ds*(0.1 + [0:nx-1]);
yy = center(2) - side(2)/2 + ds*(0.1 + [0:ny-1]);

[XX YY] = meshgrid(xx , yy);
ZZ = XX + i*YY;

%  Evaluate the analytic function on the grid
FZ = fun(ZZ);

Abf=abs(FZ);
asort=sort(Abf(:));  nz=length(asort); 
crop=0.05; Ncon=15;

%  Truncate the peaks 
top=asort(round(nz*(1-crop)));    bot=asort(1);
k=find(Abf > top);
Abf(k)=top*ones(length(k), 1);


%  Begin plotting by layers
figure(1)
clf

%  Bottom layer - color contours of |f(z)|
colormap(cool)
pcolor(xx, yy, Abf)
shading flat
colorbar

%  Assign contour levels
deltz = (top - bot) / Ncon; lg=log10(deltz);
n10=floor(lg); f10=lg - n10;
nice=[1 2 2.5 5 10];
[t, i10]=min(abs(f10 - log10(nice)));
deltz = nice(i10) * 10^n10;
zlev=deltz*(floor(bot/deltz) + [1 : Ncon]);

%  Layer 2 - Superpose black contours of magnitude
hold on
C=contour(xx,yy, Abf, zlev, 'k');  
%  clabel(C);


% Layer 3 - Add phase contours
contour(xx,yy, imag(FZ), [0 0], 'w--');
for theta=30 : 30 : 150
  rot=exp(i*theta*pi/180);
  [C,H]=contour(xx,yy, theta+imag(rot*FZ), [theta theta], 'w');
% if (theta == 0) clabel(C,H); end
end

title ...
('Magnitude (black) & phase (white) contours of f(z) in complex z plane')
xlabel('Re z')
ylabel('Im z')
axis equal
axis([xx(1) xx(nx) yy(1) yy(ny)])

fprintf('Magnitude contours run from %g to %g in steps of %g\n', ...
         zlev(1), zlev(Ncon), deltz)
fprintf('Phase contours every 30 degrees.  Dashed=0 or 180.\n\n');


%  Layer 4 - Display  branch cuts and singular points
%  Normalized approximate Laplacian
Lap=abs(del2(real(FZ))) + abs(del2(imag(FZ))) ./ (0.00001+abs(FZ));

J = find(Lap > 0.02);
[ii jj]=ind2sub(size(Lap), J);
xo=xx(1) + ds*(jj-1);
yo=yy(1) + ds*(ii-1);

%  Examine singularities of f(z) is some detail.
tx = ds*[-5:1:5]/7;
[XX, YY]=meshgrid(tx, tx);
nk=length(xo);
for k = 1: nk
  ZZ = xo(k)+XX + i*(yo(k)+YY);
  fz = fun(ZZ);
  lap=abs(del2(real(fz))) + abs(del2(imag(fz))) ./ (0.000001+abs(fz));
  J = find(lap > .105);
  [ii jj]=ind2sub(size(lap), J);
  xq=xo(k) + tx(1) + (jj-1)*(tx(2)-tx(1));
  yq=yo(k) + tx(1) + (ii-1)*(tx(2)-tx(1));
  plot(xq,yq,'y.', xq,yq,'r')
end
hold off


%  Clean up some of the large arrays
clear asort Lap


return


%  %  A small collection of analytic functions, to be put
%  %  in a separate file for execution by 'zmap'
%  function f = fun(z)
%  
%  f=exp(-sqrt(1+z.^2)) ./ (1 + 0.4*z.^2);
%
%
%  %f=(1 + z.^2).^(-1);
%  %f=exp(i*z) .* cosh(z).^(-1);
%  %f=exp(z.^(-1));
%  
%  
%  %f=tan(z);
%  %f=sqrt(z.^2 + exp(i*z));
%  %f=exp(i*z) .* cosh(z).^(-1);
%  
%  %f=(i - z).^(-1);
%  %f=log(z);
%  %f=sqrt(z).^(-1);
%  %f=sqrt(z.^2 + exp(i*z));
%  
%  
%  %f=z .^ 0.1;
%  %f=tan(z);
%  %f=cosh(z) ./ cosh( sinh(z));


%---------------------------------------------X here
function y = fun(z)



%=================================================
% tan
y=tan(z); return

y= z.^(-1) .* besselj(1,z) .* besselh(1,z);
return


%  Hankel function of order 0
y = besselh(0, z); return


% Simple Fourier transform
y = exp(i*z) ./(1+z.^2); return

%  Hankel function of order 0
y = besselh(0, -z); return


%  A simple pole
y=1 ./ z;
return

%  Inverse tan
y = atan(z); return

%  The log
y = log(z); return

%  A simple pole
y=1 ./ z;
return


%  The log
y = log(z); return
% Fourier transform with cut
y = exp(i*z) ./ (z .* sqrt(1 + z.^(-2))); return

y = 1 ./sqrt(sin(z).^4 + cos(z).^4); return
% A cubic polynomial
y = z.^3 + 4*z - 1; return

y=log(z + 0.5*i).*exp(-z.^2/2);
return

y= z.^(-1) .* besselj(1,z) .* besselh(1,z);
return
y= z.^(-1) .* besselj(1,z) .^2 
return




k=0.4;
[sn, cn,dn]=ellipz(z,k);
y=log(sn);
return

%  List of analytic functions for zmap.

% A cubic polynomial
y = z.^3 + 4*z - 1; return

% A simple pole
y = 1 ./z; return

% 3 poles of order 1
y =1.0 ./( z.^3 + 4*z - 1) ; return

% The square root
y =sqrt(z); return

% A pair of branch points
y = 1.0 ./sqrt(z.^2 - 1); return

% A finite branch cut
y = -i ./(z.*sqrt(1-1.0 ./z.^2)); return

% The exponential
y = exp(-i*z); return

%  The log
y = log(z); return

%  Integrands of 5 equiavlent integrals
y = 1 ./sqrt(sin(z).^4 + cos(z).^4); return
y = 3*z.^2 ./sqrt(2*cosh(z.^3)); return
y = sqrt(2.0 ./cosh(2*z));        return
y = cosh(z) ./sqrt(1+sinh(z).^4); return
y = 1.0 ./sqrt(1+z.^4); return

%  Integrands for 3 integrals
y= (1+3*z.^2)./cosh(z+z.^3); return
y = 1.0 ./cosh(z); return
y = 1.0 ./(1+z.^2);  return

%  Hankel function of order 1
y = besselh(1, -z); return

%  Integrand for a Fourier transform
y=exp(-i*z) ./sqrt(1-z.^2); return

%  Klein's invariant for the cube
y = (z.^8+14*z.^4+1).^3 ./(z.^4.*(z.^4-1).^4)/108;
return

%  Keyhole contour
y = log(1-z) ./(2*pi*(1+z.^3)); return

% A Hankel function
y = besselh(0, z); return

% Branch cut and pole
y =1.0 ./ (-0.2   -1*i +  sqrt(z-1) .* sqrt(z+1)); return


%Gaussian
y = exp(-z.^2); return

%  log(1+z)
y=log(1+z); return

% Another popular Fourier transform
y = exp(i*z) ./ cosh(z); return

% log and some approximations
y=z .*(1 - z.*(0.5 - z.*(0.333333 - z/4))); return
y=log(1+z); return
%  Continued fraction approx for ln(1+z)
y=z ./ (1+z ./(2+z ./(3+4*z ./(4+ 4*z/5)))); return

% Walter Munk's integrand
k=z; r=1; t=2; beta=0.5;
F=k.^3 ./ ((k.^2 + 1).*((beta*k).^2 + 1));
y=F .* bessel(0.0, r*k) .*exp(i*t*sqrt(k.^2 + 1));
return

% tan, its Taylor series and CF approx
y=tan(z); return
y=z.*(1+z.*z.*(0.333333+0.133333*z.*z)); return
y=z.*(1-0.06667*z.*z)./(1-0.4*z.*z); return



% Branch cuts and poles
y =1.0 ./ ((1+i)*(z-1) +  sqrt(z-1) .* sqrt(z+1)); return
y =1.0 ./ ((1+i)*(z-1) - i* sqrt(1-z) .* sqrt(z+1)); return
y =1.0 ./ ((1+i)*(z-1) +  sqrt(z-1) .* sqrt(z+1)); return
y =1.0 ./ (-0.2-i +  sqrt(z-1) .* sqrt(z+1)); return

%  Infinitely many branch cuts
y=sqrt(z.^2 + exp(i*z)); return


%  Aproximate Gamma function
t=z+3;
y=exp(-t + (t-0.5).*log(t)) .*(1 + 0.083333 ./t + 0.003472 ./t.^2);
y=y ./ (z.*(z+1).*(z+2)); return


% Fourier transform of 1/sqrt(z)
y = exp(-i*z) ./ sqrt(z); return


% Fourier transform with cut
y = exp(i*z) ./ (z .* sqrt(1 + z.^(-2))); return

%  S-C transform
y=z.^0.25 .* (z-1.34).^(-.25)  .* (z-3.64).^(-0.25) .* ...
  (z-4.336).^0.25; return

%  The Gaussian
y=exp(-z.^2 /2); return

% Essential singularity
y=exp(z.^(-1)); return



%---------------------------------------------X here
%  solid.m
%
%  Meshplot of |f(z)| from function in fun.m
%
%  Execute only after zmap has been run first.
%  Plots in figure 2 

figure(2)
mesh(xx, yy, Abf)
colormap cool   %  Use same colormap as the contour map
axis tight
title('Altitiude = |f(z)|')
xlabel('x = Re z')
ylabel('y = Im z')


[edge, iqot]= ...
max([sum(Abf(:,1)) sum(Abf(1,:)) sum(Abf(:,nx)) sum(Abf(ny,:))]);
azo=90*(iqot-1)+15;
view(azo, 15);


%---------------------------------------------X here
%  zview,m
%
%  Views of the real and inmaginary parts of f(z)

%  Determine which edge has the largest value, set viewing angles
e4 = [max(FZ(1,:)), max(FZ(:,nside)), max(FZ(nside,:)), max(FZ(:,1))];
[t, i4] = max(e4);
azim=90*i4 + 30;


% Truncate and plot the real part
FR=real(FZ);
K = find(FR > top);
FR(K)=top*ones(length(K),1);
K = find(FR < -top);
FR(K)=-top*ones(length(K),1);

figure(3)
clf
subplot(2,1,1)
mesh(xx,yy, FR)
title('Real part of f')
xlabel('Re z')
colormap('cool')
view([azim, 20]);

% Truncate and plot the imaginary part
FI=imag(FZ);
K = find(FI > top);
FI(K)=top*ones(length(K),1);
K = find(FI < -top);
FI(K)=-top*ones(length(K),1);

subplot(2,1,2)
mesh(xx,yy, FI)
title('Imaginary part of f')
xlabel('Re z')
colormap('cool')
view([azim, 20]);


% To get a hardcopy execute 
%            print -dpsc fig.ps
% But first pull down the file menu in the Figure 3 window and
% click on Paper Position; then adjust the size of the figure
% in the cartoon to be something sensible.
% Be aware this takes a while to complete because of Matlab's
% inefficient translator to PostScript


%---------------------------------------------X here
%  integrate.m
%
%  Perform a complex line integral:
%     I = integral on path(G) f(z) dz
%  where the analytic complex f(z) is given in the m file fun.m
%  You must invoke 'zmap' before running this procedure.
%
%  Then provide the integration path G in Figure 1 by mouse.
%  G is rendered smooth by interpolation.
%
%  G is a closed circuit if the scalar pathG = 0, the default, or is
%  open if pathG is set to any other value.  If pathG is not set, the
%  program guesses, based on the proximity of 1st and last points.
%
%  If there are only three points in a closed path, it will be
%  interpolated into a circle through those points.

%  User selects a smooth path between sample coordinates.

disp('With mouse click on points in the path of integration.');
disp('Double click last point to finish.');

%  Apply matlab mouse acquisition routine
colormap(pink)
[xz, yz]=getline(1, 'closed');
colormap(cool)

nz=length(xz);
z=xz + i*yz;
fprintf('          z                 f(z)           ')
[z , fun(z)]

%  If points repeat, remove the redundant ones
k=0;
dzmax=0;
for j = 2 : nz
  if (j<nz) dzmax=max([dzmax, abs(z(j)-z(j-1))]); end
  if (abs(z(j) - z(j-1)) < 0.01*dzmax) k=k+1; end
  xz(j-k)=real(z(j)); yz(j-k)=imag(z(j));
end
nz=nz - k;
xz=xz(1:nz); yz=yz(1:nz);


if (k > 0)
z=xz + i*yz
disp('Repeated points have been eliminated'); end

%  If pathG is not set, guess
pG=0;
if (exist('pathG')==0)
  if (abs(z(1)-z(nz-1)) > 1.5*dzmax) pG=1; end
else 
  pG=pathG;
end

if (pG == 0)
  disp('Integral performed on a closed path')
else
  nz=nz-1;
  fprintf('Integral performed from %g + %gi to %g + %gi\n',...
           xz(1),yz(1),xz(nz),yz(nz));
end

% If there are only three points on a closed path, create a circle
if (nz == 4 & pG == 0)
  z=z(1:3);
  A=[2*real(z), 2*imag(z), ones(3,1)];
  circle=A \ abs(z).^2;
  circle(3)=sqrt(circle(3) + circle(1)^2 + circle(2)^2);
  theta=2*pi*[0 : .05 : 1]';
  xz=circle(3)*cos(theta) + circle(1);
  yz=circle(3)*sin(theta) + circle(2);
  nz=length(xz);
  xz(nz)=xz(1);  yz(nz)=yz(1);
  disp('Only 3 sample points: interpolate a circular path')
end


%  There will be ns samples along the integration path:
ns=500;

%  Check: Is path is closed?
if (pG == 0)

%  Path closed - make curve periodic and match derivatives
  sz=sqrt(diff(xz).^2 + diff(yz).^2);
  sz=[0; cumsum(sz)];

  dxav=((xz(2)-xz(1))/(sz(2)-sz(1)) + ...
        (xz(nz)-xz(nz-1))/(sz(nz)-sz(nz-1)))/2;
  dyav=((yz(2)-yz(1))/(sz(2)-sz(1)) + ...
        (yz(nz)-yz(nz-1))/(sz(nz)-sz(nz-1)))/2;
  ds=sz(nz)/ns;
  ss=[sz(1) : ds : sz(nz)];
  
  sz=[-ds ; sz ; sz(nz)+ds];
  xz=[xz(1)-ds*dxav ; xz ; xz(nz)+ds*dxav];
  yz=[yz(1)-ds*dyav ; yz ; yz(nz)+ds*dyav];

%  Path is not closed
else

  sz=sqrt(diff(xz).^2 + diff(yz).^2);
  sz=[0; cumsum(sz)];

  ds=sz(nz)/ns;
  ss=[sz(1) : ds : sz(nz)];
  
end

xss=spline(sz, xz,  ss);
yss=spline(sz, yz,  ss);

%  Plot the smoothed path of integration on same scale as z-map
 boxx=center(1) + 0.5*side(1)* [-1,1,1,-1,-1]';
 boxy=center(2) + 0.5*side(2)* [-1,-1,1,1,-1]';

figure(2)
plot(xz, yz,'+',  xss,yss,'red', boxx,boxy,'black')
if (pG == 0)
  hold on; plot(xz(1),yz(1), 'o'); hold off
end
title('Smooth integration path (red)');
axis equal
disp('...Return cursor to Matlab text screen ...');
disp('Then press any key to continue integration.');
pause


%  Evaluate the function on the path G
ff=fun(xss + i*yss);
nf=length(xss);
fdz=(ff(2:nf)+ ff(1:nf-1)).* (diff(xss) + i*diff(yss))/ds/2;
%  Plot values 
plot(ss(1:nf-1), real(fdz),'red',   ss(1:nf-1),imag(fdz),'green')
title('Value of f(z) dz/ds on path: Real red, Imag green')
xlabel('Distance on path, s')
figure(1)


%  Perform trapezoidal rule integration
Integral = sum((ff(2:nf)+ ff(1:nf-1)).* (diff(xss) + i*diff(yss)))/2




%---------------------------------------------X here
%  turn.m
%
%  Rotates the viewing angle of an object displayed in current
%  figure; rotates about a vertical axis in 30 deg steps

dt=1.0;  % Pause time between frames in sec
% Initial azimuth, deg
[edge, iqot]= ...
max([sum(Abf(:,1)) sum(Abf(1,:)) sum(Abf(:,nx)) sum(Abf(ny,:))]);
azo=90*(iqot-1)+15;


pause(1)  %  Gives time to move cursor into figure

% Rotate about vertical axis
elev=30;
for az = [0 : 30 : 360] + azo
% disp([az, elev])
  view(az, elev)
  axis vis3d
  pause (dt)
end



%---------------------------------------------X here
%  zooom.m
%
%  Move the center of Figure 1 drawn by  zmap  and change the scale.
%  The center is the current center of the map.
%
%  The scalar variable  magnify  defines the factor; default = 1.5; 
%
%  To restore defaults in  zmap , enter "clear center"
%
%  Matlab's internal  zoom  works but does not recaculate the
%  f(z) at finer detail, or allow movement off the current region.

if (exist('magnify') == 0)
   magnify=1.5;
end
fprintf('\nMagnification factor,  magnify = %g\n',magnify)

side=side / magnify;
fprintf('New window size = %g %g\n\n', side)

zmap



%---------------------------------------------X here
%  pullback.m
%
%  Move the center of Figure 1 drawn by  zmap  and change the scale.
%  The center is the current center of the map.
%  The scalar variable  magnify  defines the factor; default = 1/1.5; 
%
%  To restore defaults in  zmap , enter "clear center"
%
%  Matlab's internal  zoom  works but does not recaculate the
%  f(z) at finer detail, or allow movement off the current region.

if (exist('magnify') == 0)
   magnify=1.5;
end

fprintf('\nScale reduction factor, 1/magnify = %g\n', 1/magnify)

side=side * magnify;
fprintf('New window size = %g %g\n\n', side)

zmap



%---------------------------------------------X here
%  pan.m
%
%  Move the center of Figure 1 drawn by  zmap  
%
%  The center is chosen with the mouse
%
%  To restore defaults in zmap enter "clear center"


disp('Place mouse at new center, then double-click')
disp('To exit press <RETURN> leaving mouse on plot')

for count = 1: 25
    [x, y]=getline(1);

    if (length(x) == 0) 
      fprintf('Exiting pan\n\n'); 
      return;
    end

    center=[x(1), y(1)];
    fprintf('\nNew center = %g %g\n\n', center)

    zmap

end


%---------------------------------------------X here
%  zmap
%  Draw a contour map of the analytic function defined in
%  the function script fun.m
%  
%  solid
%  In figure 2, draw a fishnet altitude plot of the function magnitude
%  |f(z)| 
%  
%  zview
%  In figure 3 draw two altitude plots, one for the real part the other
%  the imaginary part of the function f.
%  
%  turn
%  Rotate in 30 degree increments the most recently drawn altitude
%  plot.
%  
%  zooom
%  Recompute and redraw the zmap magnitude-phase plot on a
%  larger scale, with center unaltered.
%  
%  pullback
%  Reverse the effect of zooom.
%  
%  pan
%  Recompute and redraw the zmap plot centered on a new point
%  that is within the current viewing window.  
%  
%  integrate
%  Perform a complex contour integration numerically over a path drawn 
%  on the zmap graph with the mouse.  
%  
%  fun
%  This is not a zlab command, but a required function which
%  zmap refers to in order to get the complex function under
%  investigation.  
%  
%  
disp('Enter help zlab to see a list of commands')