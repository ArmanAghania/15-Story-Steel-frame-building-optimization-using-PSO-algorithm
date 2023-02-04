
w = wmax - (wmax - wmin)*(iter/Maxiter);
VpsoC = w*VpsoC + c1*rand(nvc,np).*(XpbestC-XpsoC) + c2*rand(nvc,np).*(XgbestC-XpsoC);
% VpsoB = w*VpsoB + c1*rand(nvb,np).*(XpbestB-XpsoB) + c2*rand(nvb,np).*(XgbestB-XpsoB);
VpsoBr = w*VpsoBr + c1*rand(nvBr,np).*(XpbestBr-XpsoBr) + c2*rand(nvBr,np).*(XgbestBr-XpsoBr);
VpsoP1 = w*VpsoP1 + c1*rand(nvp1,np).*(XpbestP1-XpsoP1) + c2*rand(nvp1,np).*(XgbestP1-XpsoP1);
VpsoP2 = w*VpsoP2 + c1*rand(nvp2,np).*(XpbestP2-XpsoP2) + c2*rand(nvp2,np).*(XgbestP2-XpsoP2);

for iii = 1:nvc
    for jjj = 1:np
        if VpsoC(iii,jjj) > vmaxC
           VpsoC(iii,jjj) = vmaxC;
        end
        if VpsoC(iii,jjj) < vminC
           VpsoC(iii,jjj) = vminC;
        end
    end
end

% for iii = 1:nvb
%     for jjj = 1:np
%         if VpsoB(iii,jjj) > vmaxB
%            VpsoB(iii,jjj) = vmaxB;
%         end
%         if VpsoB(iii,jjj) < vminB
%            VpsoB(iii,jjj) = vminB;
%         end
%     end
% end

for iii = 1:nvBr
    for jjj = 1:np
        if VpsoBr(iii,jjj) > vmaxBr
           VpsoBr(iii,jjj) = vmaxBr;
        end
        if VpsoBr(iii,jjj) < vminBr
           VpsoBr(iii,jjj) = vminBr;
        end
    end
end

for iii = 1:nvp1
    for jjj = 1:np
        if VpsoP1(iii,jjj) > vmaxP1
           VpsoP1(iii,jjj) = vmaxP1;
        end
        if VpsoP1(iii,jjj) < vminP1
           VpsoP1(iii,jjj) = vminP1;
        end
    end
end

for iii = 1:nvp2
    for jjj = 1:np
        if VpsoP2(iii,jjj) > vmaxP2
           VpsoP2(iii,jjj) = vmaxP2;
        end
        if VpsoP2(iii,jjj) < vminP2
           VpsoP2(iii,jjj) = vminP2;
        end
    end
end
