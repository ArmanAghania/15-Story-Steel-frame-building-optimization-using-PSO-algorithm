
A   = gpuArray(0.35);
Cd  = gpuArray(1.25);
H   = gpuArray(NS*NB);
Tem = gpuArray(Cd*0.05 * H ^ 0.75);
Ts  = gpuArray(0.7);
T0  = gpuArray(0.15);
S   = gpuArray(1.75);
S0  = gpuArray(1.1);
I   = gpuArray(1.2);

T = gpuArray(Te);

for iii = 1: length(T)
    if T(iii) <= Ts
        N(iii) = 1;
    elseif Ts < T(iii) && T(iii) <= 4
        N(iii) = (0.7/(4-Ts))*(T(iii) - Ts) + 1;
    elseif T(iii) > 4
        N(iii) = 1.7;
    end

    if T(iii) <= T0
        B1(iii) = S0 + (S - S0 + 1)*(T(iii)/T0);
    elseif T(iii) > T0 && T(iii) <= Ts
        B1(iii) = S + 1;
    elseif T(iii) > Ts
        B1(iii) = (S + 1)*(Ts/T(iii));
    end
end

BB = gpuArray(B1.*N);
Sa = gpuArray(A * BB * I);
Sa5050 = gpuArray(((72.134/475)^0.44)*Sa);
Sa1050 = gpuArray(Sa);
Sa250  = gpuArray(1.5*Sa);

C0 = gpuArray(1.4);
C1 = gpuArray(1.0);
C2 = gpuArray(1.0);

dt5050 = gpuArray(C0*C1*C2*Sa5050 * (Te^2)/(4*pi^2) * 9.80665019982);
dt1050 = gpuArray(C0*C1*C2*Sa1050 * (Te^2)/(4*pi^2) * 9.80665019982);
dt250  = gpuArray(C0*C1*C2*Sa250  * (Te^2)/(4*pi^2) * 9.80665019982);