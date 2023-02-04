%Modal Analysis

fid = fopen('ModalSCBF53.txt','wt');
fprintf(fid,'wipe \r\n');
fprintf(fid,'wipeAnalysis \r\n');
fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3 \r\n');
fprintf(fid,' \r\n');

fprintf(fid,'source HSSsection.tcl \r\n');
fprintf(fid,'source ISection.tcl \r\n');
fprintf(fid,' \r\n');
% 
fprintf(fid,'set DL 2500 \r\n');
fprintf(fid,'set LL 1000 \r\n');
fprintf(fid,'set m [expr ($DL + 0.25*$LL)*%d/9.80665019982] \r\n', Lb);
fprintf(fid,' \r\n');


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

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'fixY 0.0 1 1 0 \r\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

for ii = 1:NS*3
    fprintf(fid,'equalDOF %d %d 1 \r\n', Dfrqm(ii,1), Dfrqm(ii,2));
    if ii/3 ==1
       Cntr = Cntr + 1;
    end
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');


fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

for ii = 5:size(MainNodes,1)
    if MainNodes(ii,2) == 0 || MainNodes(ii,2) == 15
        fprintf(fid,'mass %d [expr $m/2] 1e-12 1e-12 \r\n',MainNodes(ii,1));
    else
        fprintf(fid,'mass %d $m 1e-12 1e-12 \r\n',MainNodes(ii,1));
    end
end

fprintf(fid,'geomTransf Linear 1 \r\n');

fprintf (fid,  'uniaxialMaterial Steel01   1   35153481   %d     0.03	 \r\n', E );


fprintf(fid,'\r\n');
fprintf(fid,'\r\n');
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
    fprintf(fid,'HSSsection %d 1 %d %d 4 4 4 4 \r\n', BrIndex, BrSection(ii,7), BrSection(ii,6));
    BrST(iii,1) = BrIndex;
    BrIndex = BrIndex + 1;
end

fprintf(fid,'set nn 10 \r\n');
%Columns

Cntr = 1;
Cntri = 1;
for iii = 1:NS
    for jjj = 1:NB + 1
        if jjj == 1 || jjj == NB +1 
            fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(Cntri,is) + 100);
            Cntr = Cntr + 1;
        else
            fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(Cntri + nvc/2,is) + 100);
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
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(1,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(1,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(4,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 3 || iii == 4
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(2,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(2,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(5,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 5
%             if      jjj == 1
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(3,is) + 100);
%             elseif  jjj == 4
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(3,is) + 100);
%             else
%                fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d  1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), XpsoC(6,is) + 100);
%             end
%             Cntr = Cntr +1;
%         end
%     end
% end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Beams

fprintf(fid,'geomTransf Linear 2 \r\n');

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB
       fprintf(fid,'element nonlinearBeamColumn %d %d %d $nn %d 2 \r\n', BE(Cntr,1), BE(Cntr,2), BE(Cntr,3), XpsoB(1,is) + 200);
       Cntr = Cntr + 1;
    end
end

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


% Cntr = 1;
% for iii = 1:NS
%     for jjj = 1:NB-1
%             if jjj == 1
%                 if     XpsoP1(iii,is) == 0
%                        continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(1,1),      BX(1,2),       BX(1,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(6,1),      BX(6,2),       BX(6,3),   XpsoBr(iii,is) + 300, m);
%                        CorE(Cntr,1) = BX(1,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(6,1);
%                        Cntr = Cntr + 1;
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(4,1),      BX(4,2),       BX(4,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(3,1),      BX(3,2),       BX(3,3),   XpsoBr(iii,is) + 300, m);
%                        CorE(Cntr,1) = BX(4,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(3,1);
%                        Cntr = Cntr + 1;
%                 else
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(1,1),      BX(1,2),       BX(1,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(3,1),      BX(3,2),       BX(3,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(4,1),      BX(4,2),       BX(4,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(6,1),      BX(6,2),       BX(6,3),   XpsoBr(iii,is) + 300, m);
%                        CorE(Cntr,1) = BX(1,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(3,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(4,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(6,1);
%                        Cntr = Cntr + 1;
%                 end
%             else
% 
%                 if     XpsoP2(iii,is) == 0
%                         continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(2,1),      BX(2,2),       BX(2,3),   XpsoBr(iii,is) + 300, m);
%                        fprintf(fid,'element corotTrussSection %d %d %d %d %d\r\n', BX(5,1),      BX(5,2),       BX(5,3),   XpsoBr(iii,is) + 300, m);
%                        CorE(Cntr,1) = BX(2,1);
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = BX(5,1);
%                        Cntr = Cntr + 1;                       
%                 end
%             end
%     end
%     BX(:,:) = BX(:,:) + 100;
%     m = m + 1;
% end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Gravity Loads

fprintf (fid,  'pattern Plain 1 Linear { \r\n');
GLoad = -(2500+0.25*1000);

for ii = 1:size(BE)
    fprintf(fid,'eleLoad -ele %d -type -beamUniform %d \r\n', BE(ii,1), GLoad);
end

fprintf (fid,  '} \n\n');

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf (fid,'constraints Plain\r\n');
fprintf (fid,'numberer Plain\r\n');
fprintf (fid,'system SparseGEN \r\n');
fprintf (fid,'test NormDispIncr   1.0e-6  60\r\n');
fprintf (fid,'algorithm Newton \r\n');
fprintf (fid,'integrator LoadControl 0.1 \r\n');
fprintf (fid,'analysis Static \r\n');
fprintf (fid,'analyze 10 \r\n');

fprintf (fid,  'loadConst  -time  0.0  \n\n');

fprintf (fid,  'set Pi [expr 4*atan(1)] \r\n');

% fprintf (fid,  'set omega2List [eigen 0]\n\n');
% fprintf (fid,  'puts  "First   Mode Period is: $omega2List"  \r\n');
% fprintf(fid,'set omega2 [lindex $omega2List 0]');
% fprintf(fid,'set NLP [expr 2*Pi/sqrt($omega2)]');
% fprintf(fid,'set i 0');
% fprintf(fid,'set NToutfile [opne NonlinearResults/NonlinearPeriods.out w]; \n\n');
% fprintf(fid,'puts $NToutfile "expr $NLP]" \r\n');
% fprintf(fid,'close $NToutfile \n\n');

% fprintf (fid,  'set eig [eigen 5]\n\n');
% fprintf (fid,  'puts  "First   Mode Period is: [expr 2*$Pi/pow([lindex $eig 0],0.5)]"  \r\n');
% fprintf (fid,  'puts  "Second Mode Period is: [expr 2*$Pi/pow([lindex $eig 1],0.5)]"  \n\n');
% 
% fprintf (fid,  'puts  "Second Mode Period is: [expr 2*$Pi/pow([lindex $eig 1],0.5)]"  \n\n');
% 
% fprintf (fid,  'set  NT1  [expr 2*$Pi/pow([lindex $eig 0],0.5)];  \r\n');
% fprintf (fid,  'set  NT2  [expr 2*$Pi/pow([lindex $eig 1],0.5)];  \r\n');
% fprintf (fid,  'set  NT3  [expr 2*$Pi/pow([lindex $eig 2],0.5)];  \r\n');
% fprintf (fid,  'set  NT4  [expr 2*$Pi/pow([lindex $eig 3],0.5)];  \r\n');
% fprintf (fid,  'set  NT5  [expr 2*$Pi/pow([lindex $eig 4],0.5)];  \n\n');
% % fprintf (fid,  'set  NT6  [expr 2*$Pi/pow([lindex $eig 5],0.5)];  \n\n');
% % fprintf (fid,  'set  NT7  [expr 2*$Pi/pow([lindex $eig 6],0.5)];  \n\n');
% % fprintf (fid,  'set  NT8  [expr 2*$Pi/pow([lindex $eig 7],0.5)];  \n\n');
% % fprintf (fid,  'set  NT9  [expr 2*$Pi/pow([lindex $eig 8],0.5)];  \n\n');
% % fprintf (fid,  'set  NT10 [expr 2*$Pi/pow([lindex $eig 9],0.5)];  \n\n');
% 
% fprintf (fid,  'set NToutfile [open NonlinearResults/NonlinearPeriods.out w];  \n\n');
% 
% fprintf (fid,  'puts $NToutfile "[expr $NT1]"  \r\n');
% fprintf (fid,  'puts $NToutfile "[expr $NT2]"  \r\n');
% fprintf (fid,  'puts $NToutfile "[expr $NT3]"  \r\n');
% fprintf (fid,  'puts $NToutfile "[expr $NT4]"  \r\n');
% fprintf (fid,  'puts $NToutfile "[expr $NT5]"  \n\n');
% % fprintf (fid,  'puts $NToutfile "[expr $NT6]"  \n\n');
% % fprintf (fid,  'puts $NToutfile "[expr $NT7]"  \n\n');
% % fprintf (fid,  'puts $NToutfile "[expr $NT8]"  \n\n');
% % fprintf (fid,  'puts $NToutfile "[expr $NT9]"  \n\n');
% % fprintf (fid,  'puts $NToutfile "[expr $NT10]" \n\n');
% 
% fprintf (fid,  'close $NToutfile  \n\n');
% 
% 
% fprintf (fid,  'recorder Node  -file  NonlinearResults/NonlinearEigen.out -node');
% Cntr = 5;
% for iii  = 1:NS
%     fprintf(fid,' %d',MainNodes(Cntr,1));
%     Cntr = Cntr + (NB +1);
% end
% fprintf(fid,' -dof 1 eigen 0 \r\n');
% 
% fprintf (fid,  'record  \n\n');
% 
% fprintf (fid,'puts "End of Modal5_nonlinear analyze" \r\n');
% fprintf(fid,'exit \r\n');
% 
% fclose (fid);
% 
% [status,~] = system('opensees ModalSCBF53.txt');
% 
% load NonlinearResults/NonlinearPeriods.out
% Te = NonlinearPeriods(1);
% 
% load NonlinearResults/NonlinearEigen.out
% First_Mode_Shape = NonlinearEigen(1,:)/NonlinearEigen(1,end);

fprintf (fid,  'set Pi [expr 4*atan(1)] \r\n');

fprintf (fid,  'set eig [eigen 5]\n\n');
fprintf (fid,  'puts  "First   Mode Period is: [expr 2*$Pi/pow([lindex $eig 0],0.5)]"  \r\n');
fprintf (fid,  'puts  "Second Mode Period is: [expr 2*$Pi/pow([lindex $eig 1],0.5)]"  \n\n');

fprintf (fid,  'set  NT1  [expr 2*$Pi/pow([lindex $eig 0],0.5)];  \r\n');
fprintf (fid,  'set  NT2  [expr 2*$Pi/pow([lindex $eig 1],0.5)];  \r\n');
fprintf (fid,  'set  NT3  [expr 2*$Pi/pow([lindex $eig 2],0.5)];  \r\n');
fprintf (fid,  'set  NT4  [expr 2*$Pi/pow([lindex $eig 3],0.5)];  \r\n');
fprintf (fid,  'set  NT5  [expr 2*$Pi/pow([lindex $eig 4],0.5)];  \n\n');

fprintf (fid,  'set NToutfile [open NonlinearResults/NonlinearPeriods.out w];  \n\n');

fprintf (fid,  'puts $NToutfile "[expr $NT1]"  \r\n');
fprintf (fid,  'puts $NToutfile "[expr $NT2]"  \r\n');
fprintf (fid,  'puts $NToutfile "[expr $NT3]"  \r\n');
fprintf (fid,  'puts $NToutfile "[expr $NT4]"  \r\n');
fprintf (fid,  'puts $NToutfile "[expr $NT5]"  \n\n');

fprintf (fid,  'close $NToutfile  \n\n');

fprintf (fid,  'recorder Node  -file  NonlinearResults/NonlinearEigen.out -node');
Cntr = 5;
for iii  = 1:NS
    fprintf(fid,' %d',MainNodes(Cntr,1));
    Cntr = Cntr + (NB +1);
end
fprintf(fid,' -dof 1 eigen 0 \r\n');

% fprintf (fid,  'recorder Node  -file  NonlinearResults/NonlinearEigen.out -node 110 210 310 410 510 610 710 810 910 1010 -dof 1 eigen 0  \r\n');
fprintf (fid,  'record  \n\n');

fprintf (fid,'puts "End of Modal5_nonlinear analyze" \r\n');
fprintf(fid,'exit \r\n');

fclose (fid);

[status,~] = system('opensees ModalSCBF53.txt');

load NonlinearResults/NonlinearPeriods.out
Te = gpuArray(NonlinearPeriods(1));

load NonlinearResults/NonlinearEigen.out
First_Mode_Shape = gpuArray(NonlinearEigen(1,:)/NonlinearEigen(1,end));
