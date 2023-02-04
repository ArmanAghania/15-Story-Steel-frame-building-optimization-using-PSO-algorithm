%Pushover Analysis

fid = fopen('PushoverSCBF53.txt','wt');
fprintf(fid,'wipe \r\n');
fprintf(fid,'wipeAnalysis \r\n');
fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3 \r\n');
fprintf(fid,' \r\n');

fprintf(fid,'source HSSsection.tcl \r\n');
fprintf(fid,'source ISection.tcl \r\n');
fprintf(fid,' \r\n');

fprintf(fid,'set DL 2500 \r\n');
fprintf(fid,'set LL 1000 \r\n');
fprintf(fid,'set m [expr ($DL + 0.25*$LL)*%d/9.80665019982] \r\n', Lb);
fprintf(fid,' \r\n');

fprintf(fid,'\r\n');

%Define Nodes

for ii = 1:size(MainNodes,1)
		fprintf(fid,'node %d %d %d  \r\n', MainNodes(ii,1), MainNodes(ii,2), MainNodes(ii,3));
end

fprintf(fid,'\r\n');

for ii = 7:size(EDNodes,1)
        fprintf(fid,'node %d %d %d \r\n', EDNodes(ii,1), EDNodes(ii,2), EDNodes(ii,3));
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Diaphragm
Cntr = 1;
Dfrqm = zeros(NS*3,2,'gpuArray');
for iii = 1:NS
    for jjj = 2:4
        Dfrqm(Cntr,1) = iii*100+10;
        Dfrqm(Cntr,2) = iii*100+(jjj*10);
        Cntr = Cntr + 1;
    end
end

for ii = 1:NS*3
    fprintf(fid,'equalDOF %d %d 1 \r\n', Dfrqm(ii,1), Dfrqm(ii,2));
    if ii/3 ==1
       Cntr = Cntr + 1;
    end
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'fixY 0.0 1 1 0 \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'geomTransf Linear 1 \r\n');
fprintf(fid,'geomTransf PDelta 2 \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf (fid,  'uniaxialMaterial Steel01   1   35153481   %d     0.03	 \r\n', E );
% np = 10;

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Section

CST = zeros(size(CSection,1),1,'gpuArray');
CIndex = 101;
for ii = 1:size(CSection)
    fprintf(fid,'ISection %d 1 %d %d %d %d 20 2 20 4\r\n', CIndex, CSection(ii,1), CSection(ii,2), CSection(ii,4), CSection(ii,3));
    CST(iii,1) = CIndex;
    CIndex = CIndex + 1;
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

BST = zeros(size(BSection,1),1,'gpuArray');
BIndex = 201;
for ii = 1:size(BSection)
    fprintf(fid,'ISection %d 1 %d %d %d %d 20 2 20 4\r\n', BIndex, BSection(ii,2), BSection(ii,3), BSection(ii,5), BSection(ii,4));
    BST(iii,1) = BIndex;
    BIndex = BIndex + 1;
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

BrST = zeros(size(BrSection,1),1,'gpuArray');
BrIndex = 301;
for ii = 1:size(BrSection)
    fprintf(fid,'HSSsection %d 1 %d %d 4 4 4 4\r\n', BrIndex, BrSection(ii,7), BrSection(ii,6));
    BrST(iii,1) = BrIndex;
    BrIndex = BrIndex + 1;
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'set nn 10 \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Columns

Cntr = 1;
Cntri = 1;
for iii = 1:NS
    for jjj = 1:NB + 1
        if jjj == 1 || jjj == NB +1 
            fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(Cntri,is) + 100);
            Cntr = Cntr + 1;
        else
            fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(Cntri + nvc/2,is) + 100);
            Cntr = Cntr + 1;
        end
    end
    if mod(iii,2) == 0
       Cntri = Cntri + 1;
    end
end

% Cntr = 1;
% for iii = 1:NS
%     for jjj = 1:NB + 1
%         if iii == 1 || iii == 2
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(1,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(1,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(4,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 3 || iii == 4
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(2,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(2,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(5,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 5
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(3,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(3,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(6,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%     end
% end

% Cntr = 1;
% for iii = 1:NS
%     for jjj = 1:NB + 1
%         if iii == 1 || iii == 2
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(1,is),1));
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(1,is),1));
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(4,is),1));
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 3 || iii == 4
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(2,is),1));
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(2,is),1));
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(5,is),1));
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 5
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(3,is),1));
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(3,is),1));
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  2 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), np, CST(XpsoC(6,is),1));
%             end
%             Cntr = Cntr +1;
%         end
%     end
% end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Beams


Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB
       fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d 1 \r\n', BE(Cntr,1), BE(Cntr,2), BE(Cntr,3), XpsoB(1,is) + 200);
       Cntr = Cntr + 1;
    end
end

% Cntr = 1;
% for iii = 1:NS
%     for jjj = 1:NB
%        fprintf(fid,'element nonlinearBeamColumn %d %d %d %d %d  1 \r\n', BE(Cntr,1), BE(Cntr,2), BE(Cntr,3), np, BST(XpsoB(iii,is),1));
%        Cntr = Cntr + 1;
%     end
% end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Hinges

fprintf (fid,  'uniaxialMaterial Elastic 2   2000000000000    \r\n');
fprintf (fid,  'uniaxialMaterial Elastic 3   2.0e-12   \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

for iii = 1:NS
    for jjj = 1:NB + 1
        if jjj == 1 
            fprintf(fid,'element zeroLength %d %d %d -mat 2 2 3 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 2);
        elseif jjj == NB + 1
            fprintf(fid,'element zeroLength %d %d %d -mat 2 2 3 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 1);
        else
            fprintf(fid,'element zeroLength %d %d %d -mat 2 2 3 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 1);
            fprintf(fid,'element zeroLength %d %d %d -mat 2 2 3 -dir 1 2 3 \r\n',iii*100+jjj*10 + 6, iii*100+jjj*10, iii*100+jjj*10 + 2);
        end
    end
end

% fprintf (fid,  'element zeroLength    1114   110  112    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1214   210  212    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1314   310  312    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1414   410  412    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1514   510  512    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1241   120  121    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1242   120  122    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2241   220  221    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2242   220  222    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3241   320  321    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3242   320  322    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4241   420  421    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4242   420  422    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5241   520  521    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5242   520  522    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1341   130  131    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1342   130  132    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2341   230  231    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2342   230  232    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3341   330  331    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3342   330  332    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4341   430  431    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4342   430  432    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5341   530  531    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5342   530  532    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1144   140  141    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1244   240  241    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1344   340  341    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1444   440  441    -mat   2  2  3   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1544   540  541    -mat   2  2  3   -dir    1  2   3    \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Braces

FeT = zeros(nvBr,1,'gpuArray');
Fcr = zeros(nvBr,1,'gpuArray');
deltacr = zeros(nvBr,1,'gpuArray');

% Fluxural Buckling
for ii = 1:size(XpsoBr,1)
    deltay = Fy/E;
    FeT(ii,1) = (pi^2*E)./((K*Lt)./(BrSection(XpsoBr(ii,is),4))).^2;
        if Fy/FeT(ii,1) <= 2.25
            Fcr(ii,1) = Fy*(0.658.^(Fy./FeT(ii,1)));
        elseif Fy/FeT(ii,1) > 2.25
            Fcr(ii,1) = 0.877*FeT(ii,1);
        end
    deltacr(ii,1) = Fcr(ii,1)/E;

end

fprintf (fid,  'set s1p %d \r\n',Fy);
fprintf (fid,  'set e1p %d \r\n',deltay);
fprintf (fid,  'set s2p %d \r\n',Fy );
fprintf (fid,  'set e2p %d \r\n',100*deltay);
fprintf (fid,  'set s3p %d \r\n',Fy);
fprintf (fid,  'set e3p %d \r\n',200*deltay);
fprintf (fid,  'set e2n %d \r\n',-5*deltay);
fprintf (fid,  'set e3n %d \r\n',-500*deltay);

Cntr = 4;
for iii = 1:nvBr
    fprintf (fid,  'uniaxialMaterial Hysteretic %d $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \r\n',Cntr, -Fcr(iii) ,-deltacr(iii,1) ,-0.2*Fcr(iii,1), -0.2*Fcr(iii,1) );
    Cntr = Cntr + 1;
end

% fprintf (fid,  'uniaxialMaterial Hysteretic 4  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \r\n', -Fcr(1) ,-deltacr(1) ,-0.2*Fcr(1), -0.2*Fcr(1) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 5  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \r\n', -Fcr(2) ,-deltacr(2) ,-0.2*Fcr(2), -0.2*Fcr(2) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 6  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \r\n', -Fcr(3) ,-deltacr(3) ,-0.2*Fcr(3), -0.2*Fcr(3) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 7  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \r\n', -Fcr(4) ,-deltacr(4) ,-0.2*Fcr(4), -0.2*Fcr(4) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 8  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(5) ,-deltacr(5) ,-0.2*Fcr(5), -0.2*Fcr(5) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 9  $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(6) ,-deltacr(6) ,-0.2*Fcr(6), -0.2*Fcr(6) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 10 $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(7) ,-deltacr(7) ,-0.2*Fcr(7), -0.2*Fcr(7) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 11 $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(8) ,-deltacr(8) ,-0.2*Fcr(8), -0.2*Fcr(8) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 12 $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(9) ,-deltacr(9) ,-0.2*Fcr(9), -0.2*Fcr(9) );
% fprintf (fid,  'uniaxialMaterial Hysteretic 13 $s1p $e1p $s2p $e2p $s3p $e3p %d %d %d $e2n %d $e3n  0.0 0.0 0.0 0.0 \n\n', -Fcr(10),-deltacr(10),-0.2*Fcr(10),-0.2*Fcr(10));

BX(1:3,:) = TE1(1:3,:);
BX(4:6,:) = TE2(1:3,:);

m = 4;

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB-1
            if jjj == 1
                if     XpsoP1(iii,is) == 0
                       continue
                elseif XpsoP1(iii,is) == 1
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(1,1),      BX(1,2),       BX(1,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(6,1),      BX(6,2),       BX(6,3),   BrSection(XpsoBr(iii,is),1), m);
                       CorE(Cntr,1) = BX(1,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(6,1);
                       Cntr = Cntr + 1;
                elseif XpsoP1(iii,is) == 2
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(4,1),      BX(4,2),       BX(4,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(3,1),      BX(3,2),       BX(3,3),   BrSection(XpsoBr(iii,is),1), m);
                       CorE(Cntr,1) = BX(4,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(3,1);
                       Cntr = Cntr + 1;
                       
                else
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(1,1),      BX(1,2),       BX(1,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(3,1),      BX(3,2),       BX(3,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(4,1),      BX(4,2),       BX(4,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(6,1),      BX(6,2),       BX(6,3),   BrSection(XpsoBr(iii,is),1), m);
                       CorE(Cntr,1) = BX(1,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(3,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(4,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(6,1);
                       Cntr = Cntr + 1;
                end
            else

                if     XpsoP2(iii,is) == 0
                        continue
                elseif XpsoP2(iii,is) == 1
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(2,1),      BX(2,2),       BX(2,3),   BrSection(XpsoBr(iii,is),1), m);
                       fprintf(fid,'element corotTruss %d %d %d %d %d\r\n', BX(5,1),      BX(5,2),       BX(5,3),   BrSection(XpsoBr(iii,is),1), m);
                       CorE(Cntr,1) = BX(2,1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(5,1);
                       Cntr = Cntr + 1;                       
                end
            end
    end
    BX(:,:) = BX(:,:) + 100;
    m = m + 1;
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

% fprintf (fid,  'recorder display SHAPE 800 10 600 600 -wipe  \r\n');
% fprintf (fid,  'vup 0 1 0  \r\n');
% fprintf (fid,  'vpn 0 0 1  \r\n');
% fprintf (fid,  'prp 2 1.5 10 50  \r\n');
% fprintf (fid,  'display 1 5 30  \r\n');

%Gravity Loads

fprintf (fid,  'pattern Plain 1 Linear { \r\n');
GLoad = -(2500+0.25*1000);

for ii = 1:size(BE)
    fprintf(fid,'eleLoad -ele %d -type -beamUniform %d \r\n', BE(ii,1), GLoad);
end

fprintf (fid,  '} \n\n');

fprintf (fid, 'constraints Plain \r\n');
fprintf (fid, 'numberer RCM \r\n');
fprintf (fid, 'system BandGeneral \r\n');
fprintf (fid, 'test NormDispIncr 1e-6 10 \r\n');
fprintf (fid, 'algorithm Newton \r\n');
fprintf (fid, 'integrator LoadControl 0.01 \r\n');
fprintf (fid, 'analysis Static \r\n');
fprintf (fid, 'analyze 100 \r\n');

fprintf (fid,  'loadConst  -time  0.0  \r\n');       
fprintf (fid,  'pattern Plain 2 Linear { \r\n');
fprintf (fid,  'load   1510  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   1410  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   1310  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   1210  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   1110  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   1010  %d   0.0   0.0   \r\n',  First_Mode_Shape(10));
fprintf (fid,  'load   910   %d   0.0   0.0   \r\n',  First_Mode_Shape(9));
fprintf (fid,  'load   810   %d   0.0   0.0   \r\n',  First_Mode_Shape(8));
fprintf (fid,  'load   710   %d   0.0   0.0   \r\n',  First_Mode_Shape(7));
fprintf (fid,  'load   610   %d   0.0   0.0   \r\n',  First_Mode_Shape(6));
fprintf (fid,  'load   510   %d   0.0   0.0   \r\n',  First_Mode_Shape(5));
fprintf (fid,  'load   410   %d   0.0   0.0   \r\n',  First_Mode_Shape(4));
fprintf (fid,  'load   310   %d   0.0   0.0   \r\n',  First_Mode_Shape(3));
fprintf (fid,  'load   210   %d   0.0   0.0   \r\n',  First_Mode_Shape(2));
fprintf (fid,  'load   110   %d   0.0   0.0   \r\n',  First_Mode_Shape(1));
fprintf (fid,  '} \r\n');

ControlNode = MainNodes(end,1);
Nsteps = (ceil(1.5*100* dt250)/100)/dt_incr;
fprintf (fid,  'set   ControlNode  %d \r\n', ControlNode);
fprintf (fid,  'set   Dmax      %d \r\n',  ceil(1.5*100* dt250)/100 );

fprintf (fid,  'recorder Node  -file  NonlinearResults/disp1540.out  -time   -node $ControlNode    -dof   1    disp \n\n');

fprintf (fid,  'recorder Node  -file  NonlinearResults/reac10.out    -time    -node   10     -dof   1    reaction  \r\n');
fprintf (fid,  'recorder Node  -file  NonlinearResults/reac20.out    -time    -node   20     -dof   1    reaction  \r\n');
fprintf (fid,  'recorder Node  -file  NonlinearResults/reac30.out    -time    -node   30     -dof   1    reaction  \r\n');
fprintf (fid,  'recorder Node  -file  NonlinearResults/reac40.out    -time    -node   40     -dof   1    reaction  \r\n');

fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift1.out     -iNode 10      -jNode 110   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift2.out     -iNode 110     -jNode 210   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift3.out     -iNode 210     -jNode 310   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift4.out     -iNode 310     -jNode 410   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift5.out     -iNode 410     -jNode 510   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift6.out     -iNode 510     -jNode 610   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift7.out     -iNode 610     -jNode 710   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift8.out     -iNode 710     -jNode 810   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift9.out     -iNode 810     -jNode 910   -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift10.out    -iNode 910     -jNode 1010  -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift11.out    -iNode 1010    -jNode 1110  -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift12.out    -iNode 1110    -jNode 1210  -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift13.out    -iNode 1210    -jNode 1310  -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift14.out    -iNode 1310    -jNode 1410  -dof   1    -perpDirn  2  \r\n');
fprintf (fid,  'recorder  Drift  -file  NonlinearResults/drift15.out    -iNode 1410    -jNode 1510  -dof   1    -perpDirn  2  \r\n');

fprintf (fid,'recorder  Element  -file  NonlinearResults/NLCForce.out  -ele');
for iii = 1:size(CE,1)
    fprintf(fid,' %d',CE(iii,1));
end
fprintf(fid,' localForce  \r\n');

fprintf (fid,'recorder  Element  -file  NonlinearResults/NLBrForce.out -ele');
for iii = 1:size(CorEE,1)
    fprintf (fid,' %d',CorEE(iii,1));
end
fprintf(fid,' globalForce  \r\n');

fprintf (fid,'recorder  Element  -file  NonlinearResults/PlasticCRotation.out  -ele');
for iii = 1:size(CE,1)
    fprintf(fid,' %d',CE(iii,1));
end
fprintf(fid,' plasticRotation  \r\n');

fprintf (fid,'recorder  Element  -file  NonlinearResults/PlasticBrStrain.out -ele');
for iii = 1:size(CorEE,1)
    fprintf (fid,' %d',CorEE(iii,1));
end
fprintf(fid,' deformation  \r\n');



fprintf (fid,'constraints Plain\r\n');
fprintf (fid,'numberer RCM\r\n');
fprintf (fid,'system SparseGEN \r\n');
fprintf (fid,'test NormDispIncr   1.0e-8   50\r\n');
fprintf (fid,'algorithm Newton \r\n');
fprintf (fid,'integrator DisplacementControl   $ControlNode  1  %d  \r\n',  dt_incr );
fprintf (fid,'analysis Static \r\n');
fprintf (fid,'analyze [expr int(%d)] \r\n',(ceil(1.5*100* dt250)/100)/dt_incr);
fprintf (fid,'loadConst -time 0.0 \r\n');

fprintf (fid,'puts "End of Pushover5 analyze" \r\n');
% fprintf(fid,'exit \r\n');

fclose (fid);

