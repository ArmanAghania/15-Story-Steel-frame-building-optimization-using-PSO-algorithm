rng('shuffle');
XpsoC  = gpuArray(xlbC  + rand(nvc,np) *(xubC-xlbC));
% XpsoB  = xlbB  + rand(nvb,np) *(xubB-xlbB);
XpsoB = gpuArray(3*ones(nvb,np));
XpsoBr = gpuArray(xlbBr + rand(nvBr,np)*(xubBr-xlbBr));
XpsoP1 = gpuArray(Plb1  + rand(nvp1,np)*(Pub1 -Plb1));
XpsoP2 = gpuArray(Plb2  + rand(nvp2,np)*(Pub2-Plb2));



for jj = 1:np
    for ii = 1:nvc
        XpsoC(ii,jj) = round(XpsoC(ii,jj));
        if XpsoC(ii,jj) < xlbC
           XpsoC(ii,jj) = xlbC;
        end
        if XpsoC(ii,jj) > xubC
           XpsoC(ii,jj) = xubC;
        end
    end
end

for jj = 1:np
    for ii = 1:nvb
        XpsoB(ii,jj) = round(XpsoB(ii,jj));
        if XpsoB(ii,jj) < xlbB
           XpsoB(ii,jj) = xlbB;
        end
        if XpsoB(ii,jj) > xubB
           XpsoB(ii,jj) = xubB;
        end
    end
end

for jj = 1:np
    for ii = 1:nvBr
        XpsoBr(ii,jj) = round(XpsoBr(ii,jj));
        if XpsoBr(ii,jj) < xlbBr
           XpsoBr(ii,jj) = xlbBr;
        end
        if XpsoBr(ii,jj) > xubBr
           XpsoBr(ii,jj) = xubBr;
        end
    end
end

for jj = 1:np
    for ii = 1:NS
        XpsoP1 = round(XpsoP1);
        if XpsoP1 < Plb1
           XpsoP1 = Plb1;
        end
        if XpsoP1 > Pub1
           XpsoP1 = Pub1;
        end
    end
end

for jj = 1:np
    for ii = 1:NS
        XpsoP2 = round(XpsoP2);
        if XpsoP2 < Plb2
           XpsoP2 = Plb2;
        end
        if XpsoP2 > Pub2
           XpsoP2 = Pub2;
        end
    end
end

Xpso = [XpsoC;XpsoB;XpsoBr;XpsoP1;XpsoP2];

vmaxC   = 0.5*(xubC-xlbC);
vmaxB   = 0.5*(xubB-xlbB);
vmaxBr  = 0.5*(xubBr-xlbBr);
vmaxP1  = 0.5*(Pub1-Plb1);
vmaxP2  = 0.5*(Pub2-Plb2);

vminC   = -vmaxC;
vminB   = -vmaxB;
vminBr  = -vmaxBr;
vminP1  = -vmaxP1;
vminP2  = -vmaxP2;

VpsoC   = vminC   + rand(nvc,np,'gpuArray')  *(vmaxC-vminC);
VpsoB   = vminB   + rand(nvb,np,'gpuArray')  *(vmaxB-vminB);
VpsoBr  = vminBr  + rand(nvBr,np,'gpuArray') *(vmaxBr-vminBr);
VpsoP1  = vminP1  + rand(NS,np,'gpuArray')*(vmaxP1-vminP1);
VpsoP2  = vminP2  + rand(NS,np,'gpuArray')*(vmaxP2-vminP1);

Vpso = [VpsoC;VpsoB;VpsoBr;VpsoP1;VpsoP2];