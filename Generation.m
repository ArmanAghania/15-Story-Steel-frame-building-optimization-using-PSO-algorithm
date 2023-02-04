MainNodes = zeros((NS+1)*(NB+1),3,'gpuArray');
EDNodes = zeros((NS+1)*(2*NB),3,'gpuArray');
CE = zeros((NS*(NB+1)),3,'gpuArray');
BE = zeros(NS*NB,3,'gpuArray');
TE1 = zeros((NS*NB),3,'gpuArray');
TE2 = zeros((NS*NB),3,'gpuArray');
% Generate Main Nodes
Cntr = 1;
for iii = 0:NS
    for jjj = 1:NB + 1
        MainNodes(Cntr,1) = iii*100 + jjj*10;
        MainNodes(Cntr,2) = (jjj-1)*Lb;
        MainNodes(Cntr,3) = iii*Ls;
        Cntr = Cntr + 1;
    end
end

% Generate EqualDOF Nodes
Cntr = 1;
for iii = 0:NS
    for jjj = 1:NB + 1
        if jjj == 1
            EDNodes(Cntr,1) = iii*100 + jjj*10 + 2;
            EDNodes(Cntr,2) = (jjj-1)*Lb;
            EDNodes(Cntr,3) = iii*Ls;
            Cntr = Cntr + 1;
        elseif jjj == NB+1
            EDNodes(Cntr,1) = iii*100 + jjj*10 + 1;
            EDNodes(Cntr,2) = (jjj-1)*Lb;
            EDNodes(Cntr,3) = iii*Ls;
            Cntr = Cntr + 1;
        else
            EDNodes(Cntr,1) = iii*100 + jjj*10 + 1;
            EDNodes(Cntr,2) = (jjj-1)*Lb;
            EDNodes(Cntr,3) = iii*Ls;
            Cntr = Cntr + 1;

            EDNodes(Cntr,1) = iii*100 + jjj*10 + 2;
            EDNodes(Cntr,2) = (jjj-1)*Lb;
            EDNodes(Cntr,3) = iii*Ls;
            Cntr = Cntr + 1;
        end
    end
end

%Generate Column Elements

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB + 1
        CE(Cntr,1) = (iii)*100 + jjj*10 +1;
        CE(Cntr,2) = (iii-1)*100 + jjj*10;
        CE(Cntr,3) = iii*100 + jjj*10;
        Cntr = Cntr +1;
    end
end

%Generate Beam Elements

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB
        BE(Cntr,1) = iii*100 + jjj*10 + 2;
        BE(Cntr,2) = iii*100 + jjj*10 + 2;
        BE(Cntr,3) = iii*100 + (jjj+1)*10 + 1;
        Cntr = Cntr + 1;
    end
end

%Generate Truss Elements

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB
        TE1(Cntr,1) = iii*100 + jjj*10 + 3;
        TE1(Cntr,2) = (iii-1)*100 + jjj*10;
        TE1(Cntr,3) = iii*100 + (jjj+1)*10;
        Cntr = Cntr + 1;
    end
end

Cntr = 1;

for iii = 1:NS
    for jjj = 1:NB
        TE2(Cntr,1) = iii*100 + jjj*10 + 4;
        TE2(Cntr,2) = iii*100 + jjj*10;
        TE2(Cntr,3) = (iii-1)*100 + (jjj+1)*10;
        Cntr = Cntr + 1;
    end
end
