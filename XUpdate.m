
XpsonewC = ceil(XpsoC + VpsoC);
% XpsonewB = ceil(XpsoB + VpsoB);
XpsonewB  = XpsoB;
XpsonewBr = ceil(XpsoBr + VpsoBr);
XpsonewP1 = ceil(XpsoP1 + VpsoP1);
XpsonewP2 = ceil(XpsoP1 + VpsoP1);

for iii = 1:nvc
    for jjj = 1:np
        if (xlbC <= XpsonewC(iii,jjj) && XpsonewC(iii,jjj) <= xubC)
            XpsoC(iii,jjj) = XpsonewC(iii,jjj);
        end
    end
end

for iii = 1:nvb
    for jjj = 1:np
        if (xlbB <= XpsonewB(iii,jjj) && XpsonewB(iii,jjj) <= xubB)
            XpsoB(iii,jjj) = XpsonewB(iii,jjj);
        end
    end
end

for iii = 1:nvBr
    for jjj = 1:np
        if (xlbBr <= XpsonewBr(iii,jjj) && XpsonewBr(iii,jjj) <= xubBr)
            XpsoBr(iii,jjj) = XpsonewBr(iii,jjj);
        end
    end
end

for iii = 1: nvp1
    for jjj = 1:np
        if (Plb1 <= XpsonewP1(iii,jjj) && XpsonewP1(iii,jjj) <= Pub1)
            XpsoP1(iii,jjj) = XpsonewP1(iii,jjj);
        end
    end
end

for iii = 1: nvp2
    for jjj = 1:np
        if (Plb2 <= XpsonewP2(iii,jjj) && XpsonewP2(iii,jjj) <= Plb2)
            XpsoP2(iii,jjj) = XpsonewP2(iii,jjj);
        end
    end
end
