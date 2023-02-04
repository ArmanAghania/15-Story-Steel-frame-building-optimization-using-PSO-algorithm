format longG


load BMoment.out;
load BrForce.out;
load CForce.out;

C = gpuArray(1);
Cb = gpuArray(1);
K  = gpuArray(1);
Lt = gpuArray(sqrt(Ls^2 + Lb^2));
LbC = gpuArray(Ls);
LbB = gpuArray(0);
LbT = gpuArray(Lt);
Kz  = gpuArray(1);
G   = gpuArray(E/2.6);
fic = gpuArray(0.9);
fib = gpuArray(0.9);
Fy = gpuArray(35153481.31);
PuC = zeros(nvc,1,'gpuArray');
MuC = zeros(nvc,1,'gpuArray');

Cntr = 4;

for iii = 1:nvc
    PuC(iii,1) = CForce(1,Cntr);
    Cntr = Cntr + 6;
end

% PuC = [CForce(4) ;CForce(10) ;CForce(16) ;CForce(22) ;CForce(28) ;CForce(34)];

Cntr = 3;

for iii = 1:nvc
    MuC(iii,1) = max(abs([CForce(1,Cntr),CForce(1,Cntr + 3)]));
    Cntr = Cntr + 6;
end
% MuC = [max(abs([CForce(3),CForce(6)])) ;max(abs([CForce(9), CForce(12)])) ;max(abs([CForce(15),CForce(18)])) ...
%       ;max(abs([CForce(21),CForce(24)])) ;max(abs([CForce(27),CForce(30)])) ;max(abs([CForce(33),CForce(36)]))];

PuB = gpuArray(0);
MuB = gpuArray(BMoment(2)^2/(-2*Load));
PuT = zeros(size(CorEE,1),1,'gpuArray');
Cntr = 1;
for iii = 1:size(CorEE,1)
    PuT(iii,1) = BrForce(1,Cntr);
    Cntr = Cntr + 1;
end

%     [-sqrt(BrForce(1)^2+BrForce(2 )^2) ;-sqrt(BrForce(7 )^2+BrForce(8 )^2) ;-sqrt(BrForce(13)^2+BrForce(14)^2);...
%           -sqrt(BrForce(19)^2+BrForce(20)^2) ;-sqrt(BrForce(25)^2+BrForce(26)^2)];
% PuT =[-sqrt(BrForce(1 )^2+BrForce(2 )^2) ;-sqrt(BrForce(7 )^2+BrForce(8 )^2) ;-sqrt(BrForce(13)^2+BrForce(14)^2);...
%       -sqrt(BrForce(19)^2+BrForce(20)^2) ;-sqrt(BrForce(25)^2+BrForce(26)^2)];
MuT = 0;

%Column

PFSC = 0;
PnC  = zeros(nvc,1,'gpuArray');
FeCF = zeros(nvc,1,'gpuArray');
FcrCF= zeros(nvc,1,'gpuArray');
PnCF = zeros(nvc,1,'gpuArray');
MpC  = zeros(nvc,1,'gpuArray');
RtsC = zeros(nvc,1,'gpuArray');
LpC  = zeros(nvc,1,'gpuArray');
LrC  = zeros(nvc,1,'gpuArray');
MnC1 = zeros(nvc,1,'gpuArray');
MnC  = zeros(nvc,1,'gpuArray');
FcrC = zeros(nvc,1,'gpuArray');
GSC  = zeros(nvc,1,'gpuArray');
for iii = 1:size(XpsoC,1)
    if PuC(iii) >= 0
       PnC(iii,1) = Fy*CSection(XpsoC(iii,is),5);
    elseif PuC(iii) < 0
        FeCF(iii,1) = (pi^2*E)./(((K*Ls)./CSection(XpsoC(iii,is),13))^2);    %Flexural Buckling

            if Fy/FeCF(iii,1) <= 2.25
                FcrCF(iii,1) = Fy*(0.658.^(Fy./FeCF(iii,1)));
            elseif Fy/FeCF(iii,1) > 2.25
                FcrCF(iii,1) = 0.877*FeCF(iii,1);
            end
%             Cw = (CSection(XpsoC(iii,is),1)).^2.*(CSection(XpsoC(iii,is),2))^3.*(CSection(XpsoC(iii,is),4))/24;
%             FeCT(iii,1) = (pi^2*E*Cw./(Kz*Ls) + G*CSection(XpsoC(iii,is),14)).*(1./(CSection(XpsoC(iii,is),6)+CSection(XpsoC(iii,is),7)));
%     
%             if (K*Ls/CSection(XpsoC(iii,is),13)) <= 4.71*sqrt(E/Fy) || Fy/FeCT(iii) <=2.25
%                 FcrCT(iii,1) = Fy*(0.658.^(Fy./FeCT(iii)));
%             elseif (K*Ls/CSection(XpsoC(iii,is),13)) > 4.71*sqrt(E/Fy) || Fy/FeCT(iii) <=2.25
%                 FcrCT(iii,1) = 0.877*FeCT(iii); 
%             end

        PnCF(iii,1) = FcrCF(iii,1).*CSection(XpsoC(iii,is),5);
%         PnCT(iii,1) = FcrCT(iii).*CSection(XpsoC(iii,is),5);

        PnC(iii,1)  = PnCF(iii,1);
    end

        MpC(iii,1) = Fy*CSection(XpsoC(iii,is),8);
        RtsC(iii,1)= sqrt(sqrt(CSection(XpsoC(iii,is),7)*CSection(XpsoC(iii,is),18))/CSection(XpsoC(iii,is),10)); %(10-2-5-9)
        LpC(iii,1) = 1.76*(CSection(XpsoC(iii,is),13))*sqrt(E/Fy);                          %(10-2-5-6)
        LrC(iii,1) = 1.95 .*RtsC(iii,1).*(E./(0.7.*Fy)).*sqrt(CSection(XpsoC(iii,is),14).*C/(CSection(XpsoC(iii,is),10).*(CSection(XpsoC(iii,is),16)))+sqrt((CSection(XpsoC(iii,is),14).*C/(CSection(XpsoC(iii,is),10).*(CSection(XpsoC(iii,is),16)))).^2 + 6.76.*(0.7.*Fy./E).^2));  %(10-2-5-7)
    
    if LbC <= LpC(iii,1)
        MnC(iii,1) = MpC(iii,1);
    elseif LbC > LpC(iii,1) && LbC < LrC(iii,1)
        MnC1(iii,1) = Cb*(MpC(iii,1)-(MpC(iii,1)-(0.7.*Fy.*CSection(XpsoC(iii,is),10))).*((LbC-LpC(iii,1))./(LrC(iii,1)-LpC(iii,1))));
        MnC(iii,1)  = min(abs(MnC1(iii,1)),MpC(iii,1));
    elseif LbC > LrC(iii,1)
        FcrC(iii,1) = ((Cb*pi.^2.*E)/((LbC./RtsC(iii,1))^2)).*sqrt(1 + 0.078.*(CSection(XpsoC(ii,is),14).*C/(CSection(XpsoC(iii,is),10).*(CSection(XpsoC(iii,is),1)-CSection(XpsoC(iii,is),4))))*((LbC./RtsC(iii,1))^2));
        MnC1(iii,1) = FcrC(iii,1).*CSection(XpsoC(iii,is),10);
        MnC(iii,1)  = min(abs(MnC1(iii,1)),MpC(iii,1));
    end

    if abs(PuC(iii,1))/(fic*PnC(iii,1)) < 0.2
        GSC(iii,1) = (0.5)*abs(PuC(iii,1))/(fic*PnC(iii,1)) + abs(MuC(iii,1))/(fib*MnC(iii,1)) - 1;
    elseif abs(PuC(iii))/(fic*PnC(iii)) >= 0.2
        GSC(iii,1) = abs(PuC(iii,1))/(fic*PnC(iii,1)) + (8/9)*abs(MuC(iii,1))/(fib*MnC(iii,1)) - 1;
    end

    PFSC = PFSC + (max(0,GSC(iii,1)))^2;

end

%Beam

MpB = gpuArray(Fy*BSection(XpsoB(1,1),7));
MnB = gpuArray(MpB);
GSB = gpuArray((abs(MuB)/(fib*MnB)) - 1);
PFSB= gpuArray((max(0,GSB))^2);

% PFSB = 0;
% 
% for iii = 1:size(XpsoB,1)
%     if PuB >= 0
%        PnB(iii,1) = Fy*BSection(XpsoB(iii,is),1);
%     elseif PuB < 0
%         FeBF(iii,1) = (pi^2*E)./((K*Lb)./BSection(XpsoB(iii,is),13)).^2;
%         
%         if (K*Lb)/BSection(XpsoB(iii,is),13) <= 4.71*sqrt(E/Fy) || Fy/FeBF(iii) <= 2.25
%            FcrBF(iii,1) = Fy*(0.658.^(Fy./FeBF(iii)));
%         elseif (K*Lb)/BSection(XpsoB(iii,is),13) > 4.71*sqrt(E/Fy) || Fy/FeBF(iii) > 2.25
%            FcrBF(iii,1) = 0.877*FeBF(iii);
%         end
%         
%         Cw = ((BSection(XpsoB(iii,is),1)^2*(BSection(XpsoB(iii,is),2)))^3)*(BSection(XpsoB(iii,is),4))/24;
%         FeBT(iii,1) = ((pi^2*E*Cw./(Kz*Lb)) + G*BSection(XpsoB(iii,is),14).*(1./(BSection(XpsoB(iii,is),6) + BSection(XpsoB(iii,is),10))));
%         
% 
%         if (K*Lb)/BSection(XpsoB(iii,is),13) <= 4.71*sqrt(E/Fy) || Fy/FeBF(iii) <= 2.25
%             FcrBT(iii,1) = Fy*(0.658.^(Fy./FeBT(iii)));
%         elseif (K*Lb)/BSection(XpsoB(iii,is),13) > 4.71*sqrt(E/Fy) || Fy/FeBF(iii) > 2.25
%             FcrBT(iii,1) = 0.877*FeBT(iii);
%         end
%         PnBF(iii,1) = FcrBF.*BSection(XpsoB(iii,is),1);
%         PnBT(iii,1) = FcrBT.*BSection(XpsoB(iii,is),1);
% 
%         PnB(iii,1)  = min(PnBF(iii),PnBT(iii));
% 
%     end
% 
%     MpB(iii,1) = Fy*BSection(XpsoB(iii,is),7);
%     MnB(iii,1) = MpB(iii);
% 
%     if abs(PuB)/(fib*PnB(iii)) < 0.2
%         GSB(iii,1) = ((0.5)*abs(PuB)./(fib*PnB(iii)) + abs(MuB)./(fib.*MnB(iii))) - 1;
%     elseif abs(PuB)/(fib*PnB(iii)) >= 0.2
%         GSB(iii,1) = abs(PuB)./(fib.*PnB(iii)) + (8/9).*(abs(MuB)./(fib.*MnB(iii))) - 1;
%     end
%     
%     PFSB = PFSB + max(0,GSB(iii))^2;
% 
% end

%Braces

FETF = zeros(nvBr,1,'gpuArray');
FcrTF= zeros(nvBr,1,'gpuArray');
PnT  = zeros(nvBr,1,'gpuArray');


PFST = 0;
Cntri = 1;
Cntrj = 1;

for iii = 1:size(XpsoBr,1)
    while floor(CorEE(Cntri,1)/100) == Cntrj
        FETF(iii,1) = (pi^2*E)./((K*Lt)./BrSection(XpsoBr(iii,is),4)).^2;

        if Fy/FETF(iii,1) <= 2.25
            FcrTF(iii,1) = Fy*(0.658.^(Fy./FETF(iii,1)));
        elseif Fy/FETF(iii,1) > 2.25
            FcrTF(iii,1) = 0.877*FETF(iii,1);
        end
        PnT(iii,1) = FcrTF(iii,1)*BrSection(XpsoBr(iii,is),1);

        if abs(PuT(Cntri,1))/(fic*PnT(iii,1)) < 0.2
           GST(Cntri,1) = (0.5)*abs(PuT(Cntri,1))./(fic*PnT(iii,1)) - 1;
        elseif abs(PuT(Cntri,1))/(fic*PnT(iii,1)) >= 0.2
           GST(Cntri,1) = abs(PuT(Cntri,1))./(fic*PnT(iii,1)) - 1;
        end

        PFST = PFST + (max(0,GST(Cntri,1)))^2;
        Cntri = Cntri + 1;
        if Cntri > size(CorEE,1)
            break;
        end
    end
    if Cntri > size(CorEE,1)
        break;
    end
    Cntrj = Cntrj + 1;
end
PFST = PFST + (max(0,GST(iii,1)))^2;

% for iii = 1:size(XpsoBr,1)
%     if PuT(iii) >= 0
%         PnT = Fy*BrSection(XpsoBr(iii,is),1);
%     elseif PuT(iii) < 0
% 
%         FETF(iii,1) = (pi^2*E)./((K*Lt)./BrSection(XpsoBr(iii,is),4))^2;
% 
%         if K*Lt/BrSection(XpsoBr(iii,is),4) <= 4.71*sqrt(E/Fy) || Fy/FETF(iii) <= 2.25
%            FcrTF(iii,1) = Fy*(0.658.^(Fy./FETF(iii)));
%         elseif K*Lt/BrSection(XpsoBr(iii,is),4) >  4.71*sqrt(E/Fy) || Fy/FETF(iii) > 2.25
%            FcrTF(iii,1) = 0.877*FETF(iii);
%         end
% 
%         PnT(iii,1) = (FcrTF(iii))*BrSection(XpsoBr(iii,is),1);
%     end
% 
%     MpT = [10 10 10 10 10];
%     MnT = MpT;
% 
%     if abs(PuT(iii))/(fic*PnT(iii)) < 0.2
%         GST(iii,1) = (0.5)*abs(PuT(iii))./(fic*PnT(iii))+abs(MuT)./(fib*MnT(iii)) - 1;
%     elseif abs(PuT(iii))/(fic*PnT(iii)) >= 0.2
%         GST(iii,1) = abs(PuT(iii))./(fic*PnT(iii)) + (8/9)*abs(MuT)./(fib*MnT(iii)) - 1;
%     end
% 
%     PFST = PFST + max(0,GST(iii));
% 
% end

PFS = RPS*(PFSC + PFSB + PFST);
