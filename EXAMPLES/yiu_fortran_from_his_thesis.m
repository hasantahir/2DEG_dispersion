%% A Computer program to calculated the FAR-zone field
%     n = refractive index
%     l = separation distance between two interfaces
%     k = wavenumber
%     mu = permeability of free space
%     z0 = position of the source
%     r = observation point in r-direction
%     omega = frequency
%     Ez = far-field in the z-direction
%     Er = far-field in the r-direction
%     abs_Ez = amplitude of the E-z
%     abs_Em = amplitued of E-r
%     Ephi = far-field in the phi direction
%     Rm = reflection coefficient
%     Re = reflection coefficient
%     Rmi = Approximate coefficient
%     Rei = Approximate coefficient
%% Refer to THESIS 
   %    Excitation and propagation of waves betweeen two planar interfaces
   %    1973
%%
clear
[fid,msg] = fopen('magnitude.mat', 'wt');
[fid1,msg1] = fopen('fields.mat', 'wt');
A = 4;
B = .3;
K = .125663;
omega = 6e6*2*pi;
R = 750;
L = 30;
mu = 4*pi*1e-7;
KL = K*L;
cost = omega*mu/pi*1e4;
R2 = R^2;
N = sqrt(complex(A,B));
siz = 3100;
for i = 1 : 4 : 5
    Z0 = 5*i;
    KZ0 = K*Z0;
    for j = 1 : siz
        LL(j) = j - 1;
        RR(j) = sqrt(R2+LL(j)^2);
        RRkk(j) = K*RR(j);
        expr(j) = complex(cos(RRkk(j)), sin(RRkk(j)));
        C(j) = LL(j)/RR(j);
        S(j) = R/RR(j);
        CKL(j) = KL*C(j);
        CC(j) = cos(CKL(j));
        SC(j) = sin(CKL(j));
        exprc(j) = complex(CC(j),SC(j));
        CZ(j) = cos(KZ0*C(j));
        P(j) = sqrt(complex(A - S(j)^2,B));
        RM(j) = exprc(j) * (complex(A,B)*C(j) - P(j))/(2*(P(j)*CC(j) - complex(-B,A)*C(j)^SC(j)));
        RMI(j) = exprc(j) * (N*C(j) - 1)/(2 * (CC(j) - 1i*N*C(j)*SC(j)));
        RMM(j) = abs(RM(j));
        RMMI(j) = abs(RMI(j));
        fprintf(fid,' %d %f %f %f %f \n' , j, RM(j) , RMI(j) , RMM(j) , RMMI(j));
        SMUL(j) = (complex(A,B)*C(j)*CC(j) - 1i*P(j)*SC(j))/(P(j)*CC(j) - complex(-B,A)*C(j)*SC(j));
        SMMUL(j) = abs(SMUL(j));
        ASM(j) = (1i*SC(j) - C(j)*CC(j)*N)/(1i*N*C(j)*SC(j) - CC(j));
        AMSM(j) = abs(ASM(j));
        if LL == 0
            continue
        end
        RE(j) = (C(j)-P(j))*exprc(j)/(2*(C(j)*CC(j) - 1i*P(j)*SC(j)));
        REI(j) = (C(j)-N)*exprc(j)/(2*(C(j)*CC(j)- 1i*N*SC(j)));
        RME(j) = abs(RE(j));
        RMEI(j) = abs(REI(j));
        RMUL(j) = (1i*C(j)*SC(j)+P(j)*CC(j))/(C(j)*CC(j)-1i*P(j)*SC(j)); %% CHECK
        RMMUL(j) = abs(RMUL(j));
        ARM(j) = (N*CC(j) - 1i*C(j)*SC(j))/(1i*N*SC(j) + C(j)*CC(j)); %% CHECK
        AMRM(j) = abs(ARM(j));
        PEC(j) = cost*expr(j)*sin(KZ0*C(j))/RR(j);
        Ephi(j) = PEC(j)*(1+RMUL(j));
        EMPHI(j) = abs(Ephi(j));
        AEphi(j) = PEC(j)*(1+ARM(j));
        AEMPHI(j) = abs(Ephi(j));
        
        EC(j) = cost*S(j)*CZ(j)*1i*expr(j)/RR(j);
        Ez(j) = EC(j)*S(j)*(1+SMUL(j));
        Er(j) = EC(j)*C(j)*(1+SMUL(j));
        EMZ(j) = abs(Ez(j));
        EMR(j) = abs(Er(j));
        AEZ(j) = -EC(j)*C(j)*(1+ASM(j));
        AER(j) = EC(j)*S(j)*(1+ASM(j));
        AMEZ(j) = abs(AEZ(j));
        AMER(j) = abs(AER(j));
        fprintf(fid1,' %f %f %f %f \n' , Ez(j) , Er(j) , AEZ(j) , AER(j));
    end
end
close all
plot(1:siz,abs(AEZ))

