clear all
c0 = 2.997924e8;
f = 5e9;
lam = c0/f;
k0 = 2*pi/lam;
tand = 0.001;
er = 4.8*(1 - 1i*tand);
d = lam;
k1 = sqrt(er)*k0;
numZEROS = floor(d*k0*sqrt(real(er) - 1)/pi) + 1;
orderOfTS = 0;
correctZEROS = 0;
% g = chebfun(@(v) er*sqrt(v.^2 - k0^2).*cos(d*sqrt( k1^2 - v.^2)) - ...  
%                sqrt(k1^2 - v.^2).*sin(d*sqrt( k1^2 - v.^2)),...
%                [-10,10]);
% % f = er*sqrt(v.^2 - k0^2).*cos(d*sqrt( k1^2 - v.^2)) - ...  
% %                sqrt(k1^2 - v.^2).*sin(d*sqrt( k1^2 - v.^2));
% rootss = roots(g,'all');
% clf;plot(imag(g))
% figure(2)
% plot(real(rootss),imag(rootss),'.')
% DTM = DTM(d,er,k0);
% syms v
% while correctZEROS < numZEROS
%     correctZEROS = 0;
%     zeroLIST = zeros(1,10);
%     if orderOfTS < numZEROS
%         orderOfTS = numZEROS;
%     else
%         orderOfTS = orderOfTS+1;
%     end
%     % 
%     Series = taylor (DTM, v,...
%         'ExpansionPoint', 1,...
%         'Order',orderOfTS);   % Create a TS polynomial of the DTM around 1
%     % Find roots 
%     testROOTS = solve(Series);  % Symbolic result
%     testROOTS = roots(testROOTS); % Numeric result
%     testROOTS = vpa(testROOTS); % Numeric conversion
%     for q = 1 : length(testROOTS)
%         if orderOfTS > 1
%             vr = testROOTS(q);
%         else
%             vr = testROOTS(1);
%         end
%         test = DTM(vr);
%         if (abs(test) <= 1 && real(vr) >= 0 ...
%                 && real(vr) <= d*k0*sqrt(real(er)-1)...
%                 && imag(vr) <= 0)
%             correctZEROS = correctZEROS + 1;
%             zeroLIST(correctZEROS) = vr;
%         end
%     end
% end
g = chebfun(@(v) er*sqrt(v.^2 - k0^2).*cos(d*sqrt( k1^2 - v.^2)) - ...  
               sqrt(k1^2 - v.^2).*sin(d*sqrt( k1^2 - v.^2)),...
               [-10,10]);
% f = er*sqrt(v.^2 - k0^2).*cos(d*sqrt( k1^2 - v.^2)) - ...  
%                sqrt(k1^2 - v.^2).*sin(d*sqrt( k1^2 - v.^2));
rootss = roots(g,'complex', 'norecursion');
clf;plot(imag(g))
figure(2)
clf;plot(real(rootss),imag(rootss),'.')
length(rootss)
            