T   = gpuArray(50);
La  = gpuArray(0.03); 
v   = gpuArray(1.0);
Mdi = [0.0 0.005 0.05 0.20 0.45 0.80 1.00]';       % MEAN DAMAGE INDEX (%)
ISD = [0.00 0.001429 0.003571 0.005 0.0009 0.015 0.02]';  % INTER STORY DRIFT RATIO LIMIT STATE (%)

T = gather(T);
drift = gather(drift);
Pbar50=(-1/T)*log(1-0.5) ;               MD_IO = max(drift(indIO , :)); 
Pbar10=(-1/T)*log(1-0.1)  ;               MD_LS = max(drift(indLS , :));
Pbar2= (-1/T)*log(1-0.02) ;               MD_CP = max(drift(indCP , :)); 

Pba=[Pbar50;Pbar10;Pbar2];       % probability of exeedance
MD=[MD_IO;MD_LS;MD_CP];        

Drift_Fit=fit([MD_IO MD_LS MD_CP]',[Pbar50 Pbar10 Pbar2]','exp1');

Drift_Fit.a;
Drift_Fit.b;

Pi_D_i = feval(Drift_Fit,ISD);
Pi_D_i(isnan(Pi_D_i))=0;
Pdi_D = (-1/v)*log(1-Pi_D_i);
for ii = 1:length(ISD)-1
    pp_D=isreal(Pdi_D(ii,:));
    if pp_D==1
        Pdi_D(ii,:)=Pdi_D(ii,:);
    else
        Pdi_D(ii,:)= sqrt(((real(Pdi_D(ii,:)))^2)+((imag(Pdi_D(ii,:)))^2));
    end
    if Pdi_D(ii) > 1
        Pdi_D(ii) = 1;
    end
end

for ii = 1:length(ISD)-1
    P_D_i(ii) = Pdi_D(ii)-Pdi_D(ii+1);
    if P_D_i(ii) > 1
        P_D_i(ii) = 1;
    end
end

P_D_i(length(ISD)) = Pdi_D(length(ISD));
LC = ((1/La)*(1-exp(-1*La*50))*sum(P_D_i'.*Mdi))*WF;



