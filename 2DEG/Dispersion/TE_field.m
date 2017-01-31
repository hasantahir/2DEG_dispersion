function TE_mode_field = TE_field(beta,index_mesh,x,k_zero)
% Determines TE optical field for all layers
%
% x - grid created in mesh_x.m
TotalMesh = length(x); % total number of mesh points
%
zz=beta*k_zero;
%
% Creation of constants at each mesh point
kappa = 0;
for n = 1:(TotalMesh)
    kappa(n)=sqrt((k_zero*index_mesh(n))^2-zz^2);
end
%
% Establish boundary conditions in first layer (substrate).
% Values of the fields U and V are numbered by index not by
% location along x-axis.
% For visualization purposes boundary conditions are set at first point.
U(1) = 1.0;
temp = imag(kappa(1));
if(temp<0), kappa(1) = - kappa(1);
end
% The above ensures that we get a field decaying in the substrate
V(1) = kappa(1);
%
for n=2:(TotalMesh)
    cc=cos( kappa(n)*(x(n)-x(n-1)) );
    ss=sin( kappa(n)*(x(n)-x(n-1)) );
    m(1,1)=cc;
    m(1,2)=-1i/kappa(n)*ss;
    m(2,1)=-1i*kappa(n)*ss;
    m(2,2)=cc;
    %
    U(n)=m(1,1)*U(n-1)+m(1,2)*V(n-1);
    V(n)=m(2,1)*U(n-1)+m(2,2)*V(n-1);
end
%
TE_mode_field = abs(U); % Finds Abs(E)
max_value = max(TE_mode_field);
h = plot(x,TE_mode_field/max_value); % plot normalized value of TE field
% adds text on x-axix and size of x label
xlabel('x (microns)','FontSize',22);
% adds text on y-axix and size of y label
ylabel('TE electric field','FontSize',22);
set(h,'LineWidth',1.5); % new thickness of plotting lines
set(gca,'FontSize',22); % new size of tick marks on both axes
pause
close all