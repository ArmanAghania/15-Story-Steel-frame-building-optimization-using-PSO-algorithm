%PBD Constraints

load  NonlinearResults/disp1540.out;
load  NonlinearResults/reac10.out;
load  NonlinearResults/reac20.out;
load  NonlinearResults/reac30.out;
load  NonlinearResults/reac40.out;

load  NonlinearResults/drift1.out;
load  NonlinearResults/drift2.out;
load  NonlinearResults/drift3.out;
load  NonlinearResults/drift4.out;
load  NonlinearResults/drift5.out;
load  NonlinearResults/drift6.out;
load  NonlinearResults/drift7.out;
load  NonlinearResults/drift8.out;
load  NonlinearResults/drift9.out;
load  NonlinearResults/drift10.out;
load  NonlinearResults/drift11.out;
load  NonlinearResults/drift12.out;
load  NonlinearResults/drift13.out;
load  NonlinearResults/drift14.out;
load  NonlinearResults/drift15.out;

load NonlinearResults/NLCForce.out
load NonlinearResults/NLBrForce.out

load  NonlinearResults/PlasticCRotation.out;
load  NonlinearResults/PlasticBrStrain.out;

indIO = ceil(dt5050/dt_incr);
indLS = ceil(dt1050/dt_incr);
indCP = ceil(dt250/dt_incr);

% cid = [XpsoC(1,is),XpsoC(4,is),XpsoC(4,is),XpsoC(1,is),XpsoC(1,is),XpsoC(4,is),XpsoC(4,is),XpsoC(1,is),XpsoC(2,is),XpsoC(5,is),XpsoC(5,is),XpsoC(2,is),XpsoC(2,is),XpsoC(5,is),XpsoC(5,is),XpsoC(2,is),XpsoC(3,is),XpsoC(6,is),XpsoC(6,is),XpsoC(3,is)];
% cid = [1,4,4,1,1,4,4,1,2,5,5,2,2,5,5,2,3,6,6,3];

Cntr = 1;
Cntri = 1;
for iii = 1:ceil(NS/2)
    for jjj = 1:2
        cid(1,Cntr)     = Cntri;
        cid(1,Cntr + 1) = Cntri + ceil(NS/2);
        cid(1,Cntr + 2) = Cntri + ceil(NS/2);
        cid(1,Cntr + 3) = Cntri;
        Cntr  = Cntr + 4;       
    end
    Cntri = Cntri + 1;
end


Cntri = 1;
Cntrj = 1;

% for iii = 1:size(XpsoBr,1)
%     while floor(CorEE(Cntri,1)/100) == Cntrj
%         brid(iii,1) = XpsoBr(Cntri,is);
%     end
% end
for iii = 1:size(CorEE,1)
    if      CorEE(iii,1) - 100 <100 && CorEE(iii,1) - 100 > 0 
       brid(iii,1) = XpsoBr(1,is);
    elseif  CorEE(iii,1) - 200 <100 && CorEE(iii,1) - 200 > 0  
       brid(iii,1) = XpsoBr(2,is);
    elseif  CorEE(iii,1) - 300 <100 && CorEE(iii,1) - 300 > 0  
       brid(iii,1) = XpsoBr(3,is);
    elseif  CorEE(iii,1) - 400 <100 && CorEE(iii,1) - 400 > 0  
       brid(iii,1) = XpsoBr(4,is);
    elseif  CorEE(iii,1) - 500 <100 && CorEE(iii,1) - 500 > 0  
       brid(iii,1) = XpsoBr(5,is);
    elseif  CorEE(iii,1) - 600 <100 && CorEE(iii,1) - 600 > 0  
       brid(iii,1) = XpsoBr(6,is);
    elseif  CorEE(iii,1) - 700 <100 && CorEE(iii,1) - 700 > 0  
       brid(iii,1) = XpsoBr(7,is);
    elseif  CorEE(iii,1) - 800 <100 && CorEE(iii,1) - 800 > 0  
       brid(iii,1) = XpsoBr(8,is);
    elseif  CorEE(iii,1) - 900 <100 && CorEE(iii,1) - 900 > 0  
       brid(iii,1) = XpsoBr(9,is);
    elseif  CorEE(iii,1) - 1000 <100 && CorEE(iii,1) - 1000 > 0  
       brid(iii,1) = XpsoBr(10,is);
    elseif  CorEE(iii,1) - 1100 <100 && CorEE(iii,1) - 1100 > 0  
       brid(iii,1) = XpsoBr(11,is);
    elseif  CorEE(iii,1) - 1200 <100 && CorEE(iii,1) - 1200 > 0  
       brid(iii,1) = XpsoBr(12,is);
    elseif  CorEE(iii,1) - 1300 <100 && CorEE(iii,1) - 1300 > 0  
       brid(iii,1) = XpsoBr(13,is);
    elseif  CorEE(iii,1) - 1400 <100 && CorEE(iii,1) - 1400 > 0  
       brid(iii,1) = XpsoBr(14,is);
    elseif  CorEE(iii,1) - 1500 <100 && CorEE(iii,1) - 1500 > 0  
       brid(iii,1) = XpsoBr(15,is);
    end 
end

Cntr = 1;
for iii = size(brid,1):-1:1
    brid1(Cntr,1) = brid(iii,1);
    Cntr = Cntr + 1;
end
%drift

if length(drift10) < Nsteps
    PFP = Pnlty;
else
        for ind = 1:NS
            drift = [drift1 , drift2, drift3, drift4, drift5, drift6, drift7, drift8, drift9, drift10, drift11, drift12, drift13, drift14, drift15];
            DAIO  = 0.005;
            DALS  = 0.015;
            DACP  = 0.020;
    
            GDIO(ind)  = drift(indIO,ind)/DAIO - 1;
            GDLS(ind)  = drift(indLS,ind)/DALS - 1;
            GDCP(ind)  = drift(indCP,ind)/DACP - 1;
        end
    
        PFPD = max(max(GDIO,0)).^2 + max(max(GDLS,0)).^2 + max(max(GDCP,0)).^2;
    
%Column
P         = zeros(size(CE,1),1,'gpuArray');
M         = zeros(size(CE,1),1,'gpuArray');
Mcl       = zeros(size(CE,1),1,'gpuArray');
tetay     = zeros(size(CE,1),1,'gpuArray');
deltaTC   = zeros(size(CE,1),1,'gpuArray');
Pcl       = zeros(size(CE,1),1,'gpuArray');
FeCFN     = zeros(size(CE,1),1,'gpuArray');
FcrCFN    = zeros(size(CE,1),1,'gpuArray');
PclCF     = zeros(size(CE,1),1,'gpuArray');
CrotIO    = zeros(size(CE,1),1,'gpuArray');
CrotLS    = zeros(size(CE,1),1,'gpuArray');
CrotCP    = zeros(size(CE,1),1,'gpuArray');
CdefIO    = zeros(size(CE,1),1,'gpuArray');
CdefLS    = zeros(size(CE,1),1,'gpuArray');
CdefCP    = zeros(size(CE,1),1,'gpuArray');
deltaTCIO = zeros(size(CE,1),1,'gpuArray');
deltaTCLS = zeros(size(CE,1),1,'gpuArray');
deltaTCCP = zeros(size(CE,1),1,'gpuArray');
GRIO      = zeros(size(CE,1),1,'gpuArray');
GRLS      = zeros(size(CE,1),1,'gpuArray');
GRCP      = zeros(size(CE,1),1,'gpuArray');
tetaIO    = zeros(size(CE,1),1,'gpuArray');
tetaLS    = zeros(size(CE,1),1,'gpuArray');
tetaCP    = zeros(size(CE,1),1,'gpuArray');

for iii = 1:size(CE,1)
    P(iii,1) = NLCForce(indCP,6*iii-2);
    M(iii,1) = max(abs(NLCForce(indCP,[6*iii-3,6*iii])));
    Mcl(iii,1) = MnC(cid(iii));
    tetay(iii,1) = (CSection(XpsoC(cid(iii),is),8)*Ry*Fy*Ls)/(6*E*CSection(XpsoC(cid(iii),is),6))*(1-P(iii,1)/(Ry*Fy*CSection(XpsoC(cid(iii),is),5)));
    deltaTC(iii,1) = Ry*Fy*Ls/E;

    PFTCIO = 0;
    PFTCLS = 0;
    PFTCCP = 0;

%     FeCFN(iii,1) = (pi^2*E)./((K*Ls)./(CSection(cid(iii),13))^2);


%     if P(iii,1) >= 0
%         Pcl(iii,1) = Fy*CSection(cid(1,iii),5);
%     elseif P(iii,1) < 0
    FeCFN(iii,1) = (pi^2*E)./((K*Ls)./(CSection(XpsoC(cid(iii),is),13))).^2;
    if Fy/FeCFN(iii,1) <= 2.25
        FcrCFN(iii,1) = Fy*(0.658.^(Fy./FeCFN(iii,1)));
    elseif Fy/FeCFN(iii,1) > 2.25
        FcrCFN(iii,1) = 0.877*FeCFN(iii,1);
    end

%                 FeCTN(iii) = ((pi^2*E*(CSection(XpsoC(iii,is),1)^2*(CSection(XpsoC(iii,is),2))^3)*CSection(XpsoC(iii,is),4)/24)./(K*Ls)) + G*CSection(XpsoC(iii,is),14).*(1./(CSection(XpsoC(iii,is),6) + CSection(XpsoC(iii,is),7)));
%     
%                 if K*Ls/CSection(XpsoC(iii,is),13) <= 4.71*sqrt(E/Fy) || Fy/FeCTN(iii) <= 2.25
%                     FcrCTN(iii) = Fy*(0.658.^(Fy./FeCTN(iii)));
%                 elseif K*Ls/CSection(XpsoC(iii,is),13) > 4.71*sqrt(E/Fy) || Fy/FeCTN(iii) > 2.25
%                     FcrCTN(iii) = 0.877*FeCTN;
%                 end

        PclCF(iii,1) = FcrCFN(iii,1).*CSection(XpsoC(cid(iii),is),5);
%                 PclCT(iii) = FcrCTN(iii).*CSection(XpsoC(iii,is),5);

        Pcl(iii,1) = PclCF(iii,1);
    
end

for iii = 1:size(CE,1)
    CrotIO(iii,1) = max(abs(PlasticCRotation(indIO,3*iii-1:3*iii)));
    CrotLS(iii,1) = max(abs(PlasticCRotation(indLS,3*iii-1:3*iii)));
    CrotCP(iii,1) = max(abs(PlasticCRotation(indCP,3*iii-1:3*iii)));

    CdefIO(iii,1) = abs(PlasticCRotation(indIO,3*iii-2));
    CdefLS(iii,1) = abs(PlasticCRotation(indLS,3*iii-2));
    CdefCP(iii,1) = abs(PlasticCRotation(indCP,3*iii-2));

    if P(iii,1)/(-Pcl(iii,1)) < 0
        deltaTCIO(iii,1) = 0.5*deltaTC(iii,1);
        deltaTCLS(iii,1) = 6*deltaTC(iii,1);
        deltaTCCP(iii,1) = 7*deltaTC(iii,1);

        GRIO(iii,1) = (CdefIO(iii,1)/deltaTCIO(iii,1) - 1);
        GRLS(iii,1) = (CdefLS(iii,1)/deltaTCLS(iii,1) - 1);
        GRCP(iii,1) = (CdefCP(iii,1)/deltaTCCP(iii,1) - 1);
    end

    if P(iii,1)/(-Pcl(iii,1)) >= 0 && P(iii,1)/(-Pcl(iii,1)) < 0.2
       tetaIO(iii,1) = tetay(iii,1);
       tetaLS(iii,1) = 9*tetay(iii,1);
       tetaCP(iii,1) = 11*tetay(iii,1);

       GRIO(iii,1) = (CrotIO(iii,1)/tetaIO(iii,1) - 1);
       GRLS(iii,1) = (CrotLS(iii,1)/tetaLS(iii,1) - 1);
       GRCP(iii,1) = (CrotCP(iii,1)/tetaCP(iii,1) - 1);

    end

    if P(iii,1)/(-Pcl(iii,1)) >= 0.2 && P(iii,1)/(-Pcl(iii,1)) <= 0.5
       tetaIO(iii,1) = 0.25*tetay(iii,1);
       tetaLS(iii,1) = 14*(1-5/3*(P(iii,1)./Pcl(iii,1)))*tetay(iii,1);
       tetaCP(iii,1) = 17*(1-5/3*(P(iii,1)./Pcl(iii,1)))*tetay(iii,1);

       GRIO(iii,1) = (CrotIO(iii,1)/tetaIO(iii,1) - 1);
       GRLS(iii,1) = (CrotLS(iii,1)/tetaLS(iii,1) - 1);
       GRCP(iii,1) = (CrotCP(iii,1)/tetaCP(iii,1) - 1);

    end

    if P(iii,1)/(-Pcl(iii,1)) > 0.5
       GRIO(iii,1) = ((P(iii,1))/(-Pcl(iii,1))) + ((M(iii,1))/(Mcl(iii,1))) - 1;
       GRLS(iii,1) = GRIO(iii,1);
       GRCP(iii,1) = GRIO(iii,1);
    end
       
    PFTCIO = PFTCIO + (max(0,GRIO(iii,1)))^2;
    PFTCLS = PFTCLS + (max(0,GRLS(iii,1)))^2;
    PFTCCP = PFTCCP + (max(0,GRCP(iii,1)))^2;

end

PFPC = (PFTCIO + PFTCLS + PFTCCP);

%Braces

PFTBRIO = 0;
PFTBRLS = 0;
PFTBRCP = 0;

BrdefIO     = zeros(size(CorEE,1),1,'gpuArray');
BrdefLS     = zeros(size(CorEE,1),1,'gpuArray');
BrdefCP     = zeros(size(CorEE,1),1,'gpuArray');
deltaTBr    = zeros(size(CorEE,1),1,'gpuArray');
PBr         = zeros(size(CorEE,1),1,'gpuArray');
deltaTBrIO  = zeros(size(CorEE,1),1,'gpuArray');
deltaTBrLS  = zeros(size(CorEE,1),1,'gpuArray');
deltaTBrCP  = zeros(size(CorEE,1),1,'gpuArray');
GdefIO      = zeros(size(CorEE,1),1,'gpuArray');
GdefLS      = zeros(size(CorEE,1),1,'gpuArray');
GdefCP      = zeros(size(CorEE,1),1,'gpuArray');
FeBrF       = zeros(size(CorEE,1),1,'gpuArray');
FcrBrF      = zeros(size(CorEE,1),1,'gpuArray');
deltaCBr    = zeros(size(CorEE,1),1,'gpuArray');
deltaCBrIO  = zeros(size(CorEE,1),1,'gpuArray');
deltaCBrLS  = zeros(size(CorEE,1),1,'gpuArray');
deltaCBrCP  = zeros(size(CorEE,1),1,'gpuArray');

for iii = 1:size(CorEE,1)

    BrdefIO(iii,1) = abs(PlasticBrStrain(indIO,iii));
    BrdefLS(iii,1) = abs(PlasticBrStrain(indLS,iii));
    BrdefCP(iii,1) = abs(PlasticBrStrain(indCP,iii));

    deltaTBr(iii,1) = Ry*Fy*Lt/E;

    if mod(CorEE(iii,1),2) == 1
       PBr(iii,1) =  sqrt(NLBrForce(indCP,3*iii-2).^2 + NLBrForce(indCP,3*iii-1).^2);
    else
       PBr(iii,1) = -sqrt(NLBrForce(indCP,3*iii-2).^2 + NLBrForce(indCP,3*iii-1).^2);
    end

    if PBr(iii,1) >= 0
       deltaTBrIO(iii,1) = 0.5*deltaTBr(iii,1);
       deltaTBrLS(iii,1) = 8  *deltaTBr(iii,1);
       deltaTBrCP(iii,1) = 11 *deltaTBr(iii,1);

       GdefIO(iii,1) = ((BrdefIO(iii,1))/deltaTBrIO(iii,1) - 1);
       GdefLS(iii,1) = ((BrdefLS(iii,1))/deltaTBrLS(iii,1) - 1); %Why is it  Different?
       GdefCP(iii,1) = ((BrdefCP(iii,1))/deltaTBrCP(iii,1) - 1);
    end

    if PBr(iii,1) < 0
       FeBrF(iii,1) = (pi^2*E)./((K*Lt)./BrSection(brid(iii,1),4)).^2;

       if Fy/FeBrF(iii,1) <= 2.25
          FcrBrF(iii,1) = Fy*(0.658.^(Fy./FeBrF(iii,1)));
       elseif  Fy/FeBrF(iii,1) > 2.25
          FcrBrF(iii,1) = 0.877*FeBrF(iii,1);
       end
       
     deltaCBr(iii,1) = Ry*FcrBrF(iii,1)*Lt/E;

       if ((K*Lt)/BrSection(brid(iii,1),4)) <= 2.1*sqrt(E/Fy)

          deltaCBrIO(iii,1) = 0.5*deltaCBr(iii,1);
          deltaCBrLS(iii,1) = 6  *deltaCBr(iii,1);
          deltaCBrCP(iii,1) = 7  *deltaCBr(iii,1);

       elseif ((K*Lt)/BrSection(brid(iii,1),4)) > 4.2*sqrt(E/Fy)
          
          deltaCBrIO(iii,1) = 0.5*deltaCBr(iii,1);
          deltaCBrLS(iii,1) = 7  *deltaCBr(iii,1);
          deltaCBrCP(iii,1) = 9  *deltaCBr(iii,1);

       else
          deltaCBrIO(iii,1) = 0.5*deltaCBr(iii,1);
          deltaCBrLS(iii,1) = interp1([2.1*sqrt(E/Fy)  4.2*sqrt(E/Fy)],[7*deltaCBr(iii,1)  6*deltaCBr(iii,1)],(K*Lt)/(BrSection(brid(iii,1),4)));
          deltaCBrCP(iii,1) = interp1([2.1*sqrt(E/Fy)  4.2*sqrt(E/Fy)],[9*deltaCBr(iii,1)  7*deltaCBr(iii,1)],(K*Lt)/(BrSection(brid(iii,1),4)));
       end

       GdefIO(iii,1) = ((BrdefIO(iii,1))/deltaCBrIO(iii,1) - 1);
       GdefLS(iii,1) = ((BrdefLS(iii,1))/deltaCBrLS(iii,1) - 1);
       GdefCP(iii,1) = ((BrdefCP(iii,1))/deltaCBrCP(iii,1) - 1);

    end
    
    PFTBRIO = PFTBRIO + (max(0,GdefIO(iii,1)))^2;
    PFTBRLS = PFTBRLS + (max(0,GdefLS(iii,1)))^2;
    PFTBRCP = PFTBRCP + (max(0,GdefCP(iii,1)))^2;

end
    
    PFPBR = (PFTBRIO + PFTBRLS + PFTBRCP);
    
    PFP = RPPF*(PFPD + PFPC + PFPBR);

end
    
% Plotit;
