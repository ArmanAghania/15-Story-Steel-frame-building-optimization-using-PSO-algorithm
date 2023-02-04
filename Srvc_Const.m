% 
% 
% % Servicablity Constraints
% PFG1 = 0;
% PFG2 = 0;
% PFG3 = 0;
% PFG4 = 0;
% PFG5 = 0;
% PFG6 = 0;
% PFG7 = 0;
% PFG8 = 0;
% PFG9 = 0;
% 
% G1 = zeros(1,np); %Geometrical Constraints
% G2 = zeros(1,np); %Geometrical Constraints
% G3 = zeros(1,np); %Geometrical Constraints
% G4 = zeros(1,np); %Geometrical Constraints
% G5 = zeros(1,np); %Geometrical Constraints
% G6 = zeros(1,np); %Geometrical Constraints
% G7 = zeros(1,np); %Geometrical Constraints
% G8 = zeros(1,np); %Geometrical Constraints
% G9 = zeros(1,np); %Geometrical Constraints
% 
% GGG = zeros(np,4);
% 
% G1(1,is) = CSection(XpsoC(2,is),4)/CSection(XpsoC(1,is),4)-1;
% G2(1,is) = CSection(XpsoC(3,is),4)/CSection(XpsoC(2,is),4)-1;
% G3(1,is) = CSection(XpsoC(5,is),4)/CSection(XpsoC(4,is),4)-1;
% G4(1,is) = CSection(XpsoC(6,is),4)/CSection(XpsoC(5,is),4)-1;
% 
% 
% PFG1 = PFG1 + (max(0,G1(is)))^2;
% PFG2 = PFG2 + (max(0,G2(is)))^2;
% PFG3 = PFG3 + (max(0,G3(is)))^2;
% PFG4 = PFG4 + (max(0,G4(is)))^2;
% 
% if PFG1 == 0
%    GGG(is,1) = 1;
% end
% 
% if PFG2 == 0
%    GGG(is,2) = 1;
% end
% 
% if PFG3 == 0
%    GGG(is,3) = 1;
% end
% 
% if PFG4 == 0
%    GGG(is,4) = 1;
% end
% 
% 
% if CSection(XpsoC(1,is),5) <= CSection(XpsoC(4,is),5)
% 	G5(1,is) = BSection(XpsoB(1,is),5)/CSection(XpsoC(1,is),5)-1;
% else
% 	G5(1,is) = BSection(XpsoB(1,is),5)/CSection(XpsoC(4,is),5)-1;
% end
% 
% if CSection(XpsoC(2,is),5) <= CSection(XpsoC(5,is),5)
% 	G6(1,is) = BSection(XpsoB(2,is),5)/CSection(XpsoC(2,is),5)-1;
% 	G7(1,is) = BSection(XpsoB(3,is),5)/CSection(XpsoC(2,is),5)-1;
% else
% 	G6(1,is) = BSection(XpsoB(2,is),5)/CSection(XpsoC(5,is),5)-1;
%     G7(1,is) = BSection(XpsoB(3,is),5)/CSection(XpsoC(5,is),5)-1;
% end
% 
% if CSection(XpsoC(3,is),5) <= CSection(XpsoC(6,is),5)
% 	G8(1,is) = BSection(XpsoB(4,is),5)/CSection(XpsoC(3,is),5)-1;
%   	G9(1,is) = BSection(XpsoB(5,is),5)/CSection(XpsoC(3,is),5)-1;
% 
% else
% 	G8(1,is) = BSection(XpsoB(4,is),5)/CSection(XpsoC(6,is),5)-1;
%     G9(1,is) = BSection(XpsoB(5,is),5)/CSection(XpsoC(3,is),5)-1;
% end
% 
% PFG5 = PFG5 + (max(0,G5(is)))^2;
% PFG6 = PFG6 + (max(0,G6(is)))^2;
% PFG7 = PFG7 + (max(0,G7(is)))^2;
% PFG8 = PFG8 + (max(0,G8(is)))^2;
% PFG9 = PFG9 + (max(0,G9(is)))^2;
% 
% PFG = RPG*(PFG1 + PFG2 + PFG3 + PFG4 + PFG5 + PFG6 + PFG7 + PFG8 + PFG9);
% 
% 

% PFG1 = 0;
% for iii = 1:3
%     G1(iii,1) = max(BSection(XpsoB(iii,is),2)/CSection(XpsoC(iii+3,is),2) - 1);
% 
%     PFG1 = PFG1 + (max(0,G1(iii)))^2;
% end

% PFG2 = 0;
% for iii = 1:4
%     if iii < 3
%         G2(iii,1) = CSection(XpsoC(iii+1,is),2)/CSection(XpsoC(iii,is),2) - 1;
%     end
%     if iii >= 3
%         G2(iii,1) = CSection(XpsoC(iii+2,is),2)/CSection(XpsoC(iii+1,is),2)-1;
%     end
% 
%     PFG2 = PFG2 + (max(0,G2(iii)))^2;
% end

% PFG3 = 0;
% for iii = 1:4
%     if iii < 3
%         G3(iii,1) = CSection(XpsoC(iii+1,is),1)/CSection(XpsoC(iii,is),1) - 1;
%     end
%     if iii >= 3
%         G3(iii,1) = CSection(XpsoC(iii+2,is),1)/CSection(XpsoC(iii+1,is),1)-1;
%     end
% 
%     PFG3 = PFG3 + (max(0,G3(iii)))^2;
% end
% 
% PFG4 = 0;
% for iii = 1:4
%     if iii < 3
%         G4(iii,1) = CSection(XpsoC(iii+1,is),3)/CSection(XpsoC(iii,is),3) - 1;
%     end
%     if iii >= 3
%         G4(iii,1) = CSection(XpsoC(iii+2,is),3)/CSection(XpsoC(iii+1,is),3)-1;
%     end
% 
%     PFG4 = PFG4 + (max(0,G4(iii)))^2;
% end




% PFG1 = 0;
% Cntr = 1;
% for iii = 1:5
%     G1(iii,1) = BSection(XpsoB(iii,is),2)/(min(CSection(XpsoC(Cntr,is),2),CSection(XpsoC(Cntr + 3,is),2))) - 1;
%     if mod(iii,2) == 0
%         Cntr = Cntr + 1;
%     end
%     PFG1 = PFG1 + (max(0,G1(iii)))^2;
% end

PFG2 = 0;
for iii = 1:(nvc-2)
    if iii < (nvc)/2 - 1
        G2(iii,1) = CSection(XpsoC(iii+1,is),2)/CSection(XpsoC(iii,is),2) - 1;
    end
    if iii >= (nvc)/2 - 1
        G2(iii,1) = CSection(XpsoC(iii+2,is),2)/CSection(XpsoC(iii+1,is),2)-1;
    end

    PFG2 = PFG2 + (max(0,G2(iii)))^2;
end

PFG3 = 0;
for iii = 1:(nvc-2)
    if iii < (nvc)/2 - 1
        G3(iii,1) = CSection(XpsoC(iii+1,is),1)/CSection(XpsoC(iii,is),1) - 1;
    end
    if iii >= (nvc)/2 - 1
        G3(iii,1) = CSection(XpsoC(iii+2,is),1)/CSection(XpsoC(iii+1,is),1)-1;
    end

    PFG3 = PFG3 + (max(0,G3(iii)))^2;
end

PFG4 = 0;
for iii = 1:(nvc-2)
    if iii < (nvc)/2 - 1
        G4(iii,1) = CSection(XpsoC(iii+1,is),3)/CSection(XpsoC(iii,is),3) - 1;
    end
    if iii >= (nvc)/2 - 1
        G4(iii,1) = CSection(XpsoC(iii+2,is),3)/CSection(XpsoC(iii+1,is),3)-1;
    end

    PFG4 = PFG4 + (max(0,G4(iii)))^2;
end




PFG = RPG*( PFG2 + PFG3 + PFG4);
