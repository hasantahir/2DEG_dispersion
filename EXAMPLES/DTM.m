function out = DTM( d, er, k0)
% returns a function handle in 'v'
% d is the thickness of the slab
% er is the relative permittivity of the slab dielectric 
% k0 is the free-space wavenumber

k1 = sqrt(er)*k0;
% out = @(v) er*v.*cos(sqrt((d*k0)^2 * (er-1) - v.^2));
% 
% out = @(v) er*sqrt(v.^2 - k0^2).*cos(d*sqrt( k1^2 - v.^2)) - ...  
%               sqrt(k1^2 - v.^2).*sin(d*sqrt( k1^2 - v.^2));

out = @(v) er*v.*cos(sqrt((d*k0)^2 * (er-1) - v.^2)) - ...  
              sqrt((d*k0)^2 * (er-1) - v.^2).*sin(sqrt((d*k0)^2 * (er-1) - v.^2));
          
% er*v*Cos[Sqrt[(d*k0) ˆ2*(er - 1) - vˆ2]] - 
%   Sqrt[(d*k0) ˆ2*(er - 1) - vˆ2]*Sin[Sqrt[(d*k0) ˆ2*(er - 1) - vˆ2]]
           
end