
NS = 5; %Number of Stairs
NB = 3; %Number of Bays
Ls = 3; %Length of Stairs
Lb = 5; %Length of Bays
Lt = sqrt(Ls^2+Lb^2);
np = 40;    %Number of Particles
nvc = 6;    %Comlumn Groups
nvb = NS;   %Beam Groups
nvBr = NS;  %Brace Groups
nvp1 = NS;
nvp2 = NS;
nv   = nvc + nvb + nvBr + nvp1 + nvp2;
dt_incr = 0.001;

E  = 20394323844;  %kgf/m^2
Fy = 35153481.31;  %kgf/m^2
Fu = 45699526;
Ry = 1.1;
Gamma = 7850;

Pnlty = 1e6;
RPG   = 1e3;
RPS   = 1e3;
RPPF  = 1e3;
incr  = 0.002;

CSection    = load('Input/CSection.txt');            %Column Sections, n, w, A, d, bf, tw, tf, bf/2tf, h, tw
BSection    = load('Input/BSection.txt');              %Beam   Sections, n, w, A, d, bf, tw, tf, bf/2tf, h, tw
BrSection   = load('Input/BrSection.txt');             %Brace  Sections, n, w, A, d, tf

XpsoC = [1;1;1;6;2;1];
XpsoB = [3;3;3;3;3];
XpsoBr= [16;14;14;12;10];
XpsoP1= [0;0;0;0;0];
XpsoP2= [1;1;1;1;1];