clear;
close;
clc;


NS = gpuArray(15); %Number of Stairs
NB = gpuArray(3); %Number of Bays
Ls = gpuArray(3); %Length of Stairs
Lb = gpuArray(5); %Length of Bays
Lt = gpuArray(sqrt(Ls^2+Lb^2));
np = gpuArray(1000);    %Number of Particles
nvc = gpuArray(16);    %Comlumn Groups
nvb = gpuArray(NS);   %Beam Groups
nvBr = gpuArray(NS);  %Brace Groups
nvp1 = gpuArray(NS);
nvp2 = gpuArray(NS);
nv   = gpuArray(nvc + nvb + nvBr + nvp1 + nvp2);
dt_incr = gpuArray(0.002);


xlbC  = gpuArray(1);
xubC  = gpuArray(24);

xlbB  = gpuArray(1);
xubB  = gpuArray(24);

xlbBr = gpuArray(1);
xubBr = gpuArray(23);

Plb1  = gpuArray(0);
Pub1  = gpuArray(3);

Plb2  = gpuArray(0);
Pub2  = gpuArray(1);

Maxiter  = gpuArray(200);
Maxcycle = gpuArray(10);

fitness  = [];

wmax = gpuArray(0.5);
wmin = gpuArray(0);

c1 = gpuArray(2);
c2 = gpuArray(2);

E  = gpuArray(20394323844);  %kgf/m^2
Fy = gpuArray(35153481.31);  %kgf/m^2
Fu = gpuArray(45699526);
Ry = gpuArray(1.1);
Gamma = gpuArray(7850);

% Pnlty = 1e6;
% RPG   = 1e3;
% RPS   = 1e3;
% RPPF  = 1e3;
Pnlty = gpuArray(1e6);
RPG   = gpuArray(1e3);
RPS   = gpuArray(1e3);
RPPF  = gpuArray(1e3);
incr  = gpuArray(0.002);
CH    = zeros(Maxiter,Maxcycle,'gpuArray');


CSection    = gpuArray(load('Input/CSection.txt'));            %Column Sections, n, w, A, d, bf, tw, tf, bf/2tf, h, tw
BSection    = gpuArray(load('Input/BSection.txt'));              %Beam   Sections, n, w, A, d, bf, tw, tf, bf/2tf, h, tw
BrSection   = gpuArray(load('Input/BrSection.txt'));             %Brace  Sections, n, w, A, d, tf

flag1 = gpuArray(0);
flag2 = gpuArray(0);
WF = gpuArray(1e6);
LC = gpuArray(1e6);
PF   = zeros(np,1,'gpuArray');
WW   = zeros(np,1,'gpuArray');
SOFW = zeros(np,1,'gpuArray');

Generation;

T = zeros(Maxiter,Maxcycle,'gpuArray');

for ic = 1:Maxcycle
    IXV;
    for iter = 1:Maxiter
        tic
        for is = 1:np
            if flag2 == 1
                dFiles;
            end
            flag2 = 0;
%             clc;
%             Iteration = is
%             if is == 1
%                 for iii = 1:NS
%                     for jjj = 1:nvc
% %                         WF = WF + Gamma * Ls * (Nb + 1) * (CSection(XpsoC(jjj,is),5))
%                  WF =     Gamma*(Ls*(4*(CSection(XpsoC(1,is),5) + CSection(XpsoC(2,is),5) + CSection(XpsoC(4,is),5) + CSection(XpsoC(5,is),5)) + 2*(CSection(XpsoC(3,is),5) + CSection(XpsoC(6,is),5))) + ...
%                             Lb*(3*(BSection(XpsoB(1,is),1) + BSection(XpsoB(2,is),1) + BSection(XpsoB(3,is),1) + BSection(XpsoB(4,is),1) + BSection(XpsoB(5,is),1))) + ...
%                              6*Lt*(BrSection(XpsoBr(1,is),1) + BrSection(XpsoBr(2,is),1) + BrSection(XpsoBr(3,is),1) + BrSection(XpsoBr(4,is),1) + BrSection(XpsoBr(5,is),1)));
%             elseif CorEE == ture
%                 for iii = 1:size(CorEE,1)
%                     WF = Gamma* (Ls*(4*(CSection(XpsoC(1,is),5) + CSection(XpsoC(2,is),5) + CSection(XpsoC(4,is),5) + CSection(XpsoC(5,is),5)) + 2*(CSection(XpsoC(3,is),5) + CSection(XpsoC(6,is),5))) + ...
%                                  Lb*(3*(BSection(XpsoB(1,is),1) + BSection(XpsoB(2,is),1) + BSection(XpsoB(3,is),1) + BSection(XpsoB(4,is),1) + BSection(XpsoB(5,is),1))));
%                     WF = WF + Lt*CorEE(iii,4);
%                 end
%             end
            Srvc_Const;
            if PFG > 0
                disp('Constructability Problem')
                PFS = Pnlty;
                PFP = Pnlty;
            elseif PFG == 0
                Input_Service;
                BrNS;
                    if iter > 1 && flag1 == 1
                        if WF >= Gb
                           continue
                        end
                    end
                [status,~] = system('opensees SCBF53.txt');
                Strng_Const;
                if PFS > 0
                    disp('Strength Problem')
                    PFP = Pnlty;

                elseif PFS == 0
                    Input_PBD;
                    [status,~] = system('opensees PushoverSCBF53.txt');
                    PBD_Const;
                    flag2 = 1;
                    LCC;
                    if PFP >0
                      disp('PBD Problem')
                    else
                        disp('Alright')
                        flag1 = 1;
%                         Plotit;
                    end
                end
            end
            PF(is,1) = PFG + PFS + PFP;
            WW(is,1)  = WF + LC;
            SOFW(is,1) = WW(is,1)*(1+PF(is,1));

        end
     PBest;
     GBest;
     VUpdate;
     XUpdate;
     minn = find(SOFW == Gb);
     [nmin mmin] = size(minn);
%      if mmin >round(0.8*np)
%         Xgbest;
%         break
%      end
     clc
     Gb
     CH(iter,ic) = Gb;
     Monitor = [iter,ic]
     T(iter,ic) = toc;
    end
    GB(ic) = Gb;
    XGB(1:size(Xpso,1),ic) = Xgbest(:,1);
    XGB(size(Xpso,1) + 1,ic) = GB(ic);
    BestS = XGB(ic,:);
    filename = [num2str(ic) '.xlsx'];
    writetable(BestS,filename,'Sheet',ic);
end

