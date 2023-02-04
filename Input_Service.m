E   = gpuArray(20394323844);
Nu  = gpuArray(0.3);
G   = gpuArray(E / (2*(1+Nu)));
Rho = gpuArray(7850);


%OpenSees Input

fid = fopen('SCBF53.txt','w');
fprintf(fid,'wipe \r\n');
fprintf(fid,'wipeAnalysis \r\n');
fprintf(fid,'model BasicBuilder -ndm 2 -ndf 3 \r\n');
fprintf(fid,' \r\n');

% fprintf(fid,'source HSSsection.tcl \r\n');
% fprintf(fid,'source ISection.tcl \r\n');
% fprintf(fid,' \r\n');
% 
% fprintf(fid,'set DL 2500 \r\n');
% fprintf(fid,'set LL 1000 \r\n');
% fprintf(fid,'set m [expr ($DL + 0.25*$LL)*%d/9.80665019982] \r\n', Lb);
% fprintf(fid,' \r\n');

%Define Nodes

for ii = 1:size(MainNodes,1)
		fprintf(fid,'node %d %d %d  \r\n', MainNodes(ii,1), MainNodes(ii,2), MainNodes(ii,3));
end

fprintf(fid,'\r\n');

for ii = 7:size(EDNodes,1)
        fprintf(fid,'node %d %d %d \r\n', EDNodes(ii,1), EDNodes(ii,2), EDNodes(ii,3));
end

fprintf(fid,'fixY 0.0 1 1 0 \r\n');

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


fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'geomTransf Linear 1 \r\n');


%Sections
% fprintf (fid,  'uniaxialMaterial Elastic 3   2000000000000    \r\n');

% CIndex = 101;
% for ii = 1:size(CSection)
% 
%     fprintf(fid,'ISection %d 3 %d %d %d %d 5 5 5 5\r\n', CIndex, CSection(ii,1), CSection(ii,2), CSection(ii,4), CSection(ii,3));
%     CIndex = CIndex + 1;
% end
% 
% BIndex = 201;
% for ii = 1:size(BSection)
%     fprintf(fid,'ISection %d 3 %d %d %d %d 5 5 5 5\r\n', BIndex, BSection(ii,1), BSection(ii,2), BSection(ii,4), BSection(ii,3));
%     BIndex = BIndex + 1;
% end
% 
% BrIndex = 301;
% for ii = 1:size(BrSection)
%     fprintf(fid,'HSSsection %d 3 %d %d 5 2 5 2\r\n', BrIndex, BrSection(ii,7), BrSection(ii,6));
%     BrIndex = BrIndex + 1;
% end

%Define Elements

%Column
Cntr = 1;
Cntri = 1;
for iii = 1:NS
    for jjj = 1:NB + 1
        if jjj == 1 || jjj == NB +1 
            fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n',CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(Cntri,is),5), E, CSection(XpsoC(Cntri,is),6));
            Cntr = Cntr + 1;
        else
            fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n',CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(Cntri + nvc/2,is),5), E, CSection(XpsoC(Cntri + nvc/2,is),6));
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
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(1,is),5), E, CSection(XpsoC(1,is),6));
%             elseif  jjj == 4
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(1,is),5), E, CSection(XpsoC(1,is),6));
%             else
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(4,is),5), E, CSection(XpsoC(4,is),6));
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 3 || iii == 4
%             if      jjj == 1
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(2,is),5), E, CSection(XpsoC(2,is),6));
%             elseif  jjj == 4
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(2,is),5), E, CSection(XpsoC(2,is),6));
%             else
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(5,is),5), E, CSection(XpsoC(5,is),6));
%             end
%             Cntr = Cntr +1;
%         end
%         if iii == 5
%             if      jjj == 1
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(3,is),5), E, CSection(XpsoC(3,is),6));
%             elseif  jjj == 4
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(3,is),5), E, CSection(XpsoC(3,is),6));
%             else
%                fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 1 \r\n', CE(Cntr,1), CE(Cntr,2), CE(Cntr,3), CSection(XpsoC(6,is),5), E, CSection(XpsoC(6,is),6));
%             end
%             Cntr = Cntr +1;
%         end
%     end
% end
 
fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

%Beam
fprintf(fid,'geomTransf Linear 2 \r\n');

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB
       fprintf(fid,'element elasticBeamColumn %d %d %d %d %d %d 2 \r\n', BE(Cntr,1), BE(Cntr,2), BE(Cntr,3), BSection(XpsoB(iii,is),1), E, BSection(XpsoB(iii,is),6));
       Cntr = Cntr + 1;
    end
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');


%Hinges

fprintf (fid,  'uniaxialMaterial Elastic 1   2000000000000    \r\n');
fprintf (fid,  'uniaxialMaterial Elastic 2   2.0e-12   \r\n');

for iii = 1:NS
    for jjj = 1:NB + 1
        if jjj == 1 
            fprintf(fid,'element zeroLength %d %d %d -mat 1 1 2 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 2);
        elseif jjj == NB + 1
            fprintf(fid,'element zeroLength %d %d %d -mat 1 1 2 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 1);
        else
            fprintf(fid,'element zeroLength %d %d %d -mat 1 1 2 -dir 1 2 3 \r\n',iii*100+jjj*10 + 5, iii*100+jjj*10, iii*100+jjj*10 + 1);
            fprintf(fid,'element zeroLength %d %d %d -mat 1 1 2 -dir 1 2 3 \r\n',iii*100+jjj*10 + 6, iii*100+jjj*10, iii*100+jjj*10 + 2);
        end
    end
end

% fprintf (fid,  'element zeroLength    1114   110  112    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1214   210  212    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1314   310  312    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1414   410  412    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1514   510  512    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1241   120  121    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1242   120  122    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2241   220  221    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2242   220  222    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3241   320  321    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3242   320  322    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4241   420  421    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4242   420  422    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5241   520  521    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5242   520  522    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1341   130  131    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1342   130  132    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2341   230  231    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    2342   230  232    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3341   330  331    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    3342   330  332    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4341   430  431    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    4342   430  432    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5341   530  531    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    5342   530  532    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1144   140  141    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1244   240  241    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1344   340  341    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1444   440  441    -mat   1  1  2   -dir    1  2   3    \r\n');
% fprintf (fid,  'element zeroLength    1544   540  541    -mat   1  1  2   -dir    1  2   3    \r\n');


fprintf(fid,'\r\n');
fprintf(fid,'\r\n');


%Braces

fprintf (fid,  'uniaxialMaterial Elastic 4        20394323844  \r\n');


% index = 1;
% index2 = 2;
% BX = zeros(2*NB,3);
% Cntr = 1;
% for jjj = 1:2
%     for iii = 1:2*NB
%         if jjj == 1
%             BX(iii,1) = TE1(Cntr,1);
%             BX(iii,2) = TE1(Cntr,2);
%             BX(iii,3) = TE1(Cntr,3);
%             Cntr = Cntr + 1;
%         else
%             Cntr = 1;
%             BX(iii,1) = TE2(Cntr,1);
%             BX(iii,2) = TE2(Cntr,2);
%             BX(iii,3) = TE2(Cntr,3);
%             Cntr = Cntr + 1;
%         end
%     end
% end

BX(1:3,:) = TE1(1:3,:);
BX(4:6,:) = TE2(1:3,:);
       

Cntr = 1;
for iii = 1:NS
    for jjj = 1:NB-1
            if jjj == 1
                if     XpsoP1(iii,is) == 0
                       continue
                elseif XpsoP1(iii,is) == 1
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(1,1),      BX(1,2),       BX(1,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(6,1),      BX(6,2),       BX(6,3),   BrSection(XpsoBr(iii,is),1));
                       CorE(Cntr,1) = BX(1,1);
                       CorE(Cntr,2) = BX(1,2);
                       CorE(Cntr,3) = BX(1,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(6,1);
                       CorE(Cntr,2) = BX(6,2);
                       CorE(Cntr,3) = BX(6,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                elseif XpsoP1(iii,is) == 2
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(4,1),      BX(4,2),       BX(4,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(3,1),      BX(3,2),       BX(3,3),   BrSection(XpsoBr(iii,is),1));
                       CorE(Cntr,1) = BX(4,1);
                       CorE(Cntr,2) = BX(4,2);
                       CorE(Cntr,3) = BX(4,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(3,1);
                       CorE(Cntr,2) = BX(3,2);
                       CorE(Cntr,3) = BX(3,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       
                else
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(1,1),      BX(1,2),       BX(1,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(3,1),      BX(3,2),       BX(3,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(4,1),      BX(4,2),       BX(4,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(6,1),      BX(6,2),       BX(6,3),   BrSection(XpsoBr(iii,is),1));
                       CorE(Cntr,1) = BX(1,1);
                       CorE(Cntr,2) = BX(1,2);
                       CorE(Cntr,3) = BX(1,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(3,1);
                       CorE(Cntr,2) = BX(3,2);
                       CorE(Cntr,3) = BX(3,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(4,1);
                       CorE(Cntr,2) = BX(4,2);
                       CorE(Cntr,3) = BX(4,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(6,1);
                       CorE(Cntr,2) = BX(6,2);
                       CorE(Cntr,3) = BX(6,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                end
            else

                if     XpsoP2(iii,is) == 0
                        continue
                elseif XpsoP2(iii,is) == 1
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(2,1),      BX(2,2),       BX(2,3),   BrSection(XpsoBr(iii,is),1));
                       fprintf(fid,'element corotTruss %d %d %d %d 4\r\n', BX(5,1),      BX(5,2),       BX(5,3),   BrSection(XpsoBr(iii,is),1));
                       CorE(Cntr,1) = BX(2,1);
                       CorE(Cntr,2) = BX(2,2);
                       CorE(Cntr,3) = BX(2,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;
                       CorE(Cntr,1) = BX(5,1);
                       CorE(Cntr,2) = BX(5,2);
                       CorE(Cntr,3) = BX(5,3);
                       CorE(Cntr,4) = BrSection(XpsoBr(iii,is),1);
                       Cntr = Cntr + 1;                       
                end
            end
    end
    BX(:,:) = BX(:,:) + 100;
end

% Cntr = 1;
% for iii = 1:NS
%     for jjj = 1:NB-1
%         if iii == 1
%             if jjj == 1
%                 if XpsoP1(iii,is) == 0
%                     continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 113 12 121 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 134 132 41 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 113;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 134;
%                        Cntr = Cntr + 1;
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection 114 112 21 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 133 32 141 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 114;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 133;
%                        Cntr = Cntr + 1;
%                 else
%                        fprintf(fid,'element corotTrussSection 113 12 121 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 134 132 41 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 114 112 21 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 133 32 141 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 113;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 134;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 114;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 133;
%                        Cntr = Cntr + 1;
%                 end
%             elseif jjj == 2
%                 if XpsoP2(iii,is) == 0
%                     continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 123 22 131 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 124 122 31 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 123;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 124;
%                        Cntr = Cntr + 1;                       
%                 end
%             end
%         elseif iii == 2
%             if jjj == 1
%                 if XpsoP1(iii,is) == 0
%                     continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 213 112 221 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 234 232 141 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 213;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 234;
%                        Cntr = Cntr + 1;                       
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection 214 212 121 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 233 132 241 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 214;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 233;
%                        Cntr = Cntr + 1;                       
%                 else
%                        fprintf(fid,'element corotTrussSection 213 112 221 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 234 232 141 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 214 212 121 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 233 132 241 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 213;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 234;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 214;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 233;
%                        Cntr = Cntr + 1;                       
%                 end
%             elseif jjj == 2
%                 if XpsoP2(iii,is) == 0
%                     continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 223 122 231 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 224 222 131 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 223;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 224;
%                        Cntr = Cntr + 1;                                              
%                 end
%             end
%         elseif iii == 3
%             if jjj == 1
%                 if XpsoP1(iii,is) == 0
%                     continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 313 212 321 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 334 332 241 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 313;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 334;
%                        Cntr = Cntr + 1;                                              
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection 314 312 221 %d 4\r\n',    BSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 333 232 341 %d 4\r\n',    BSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 314;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 333;
%                        Cntr = Cntr + 1;                                              
%                 else
%                        fprintf(fid,'element corotTrussSection 313 212 321 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 334 332 241 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 314 312 221 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 333 232 341 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 313;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 334;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 314;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 333;
%                        Cntr = Cntr + 1;                                              
%                 end
%             elseif jjj == 2
%                 if XpsoP2(iii,is) == 0
%                     continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 323 222 331 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 324 322 231 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 323;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 324;
%                        Cntr = Cntr + 1;                                                                     
%                 end
%             end
%         elseif iii == 4
%             if jjj == 1
%                 if XpsoP1(iii,is) == 0
%                     continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 413 312 421 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 434 432 341 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 413;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 434;
%                        Cntr = Cntr + 1;                                              
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection 414 412 321 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 433 332 441 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 414;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 433;
%                        Cntr = Cntr + 1;                                                                     
%                 else
%                        fprintf(fid,'element corotTrussSection 413 312 421 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 434 432 341 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 414 412 321 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 433 332 441 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 413;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 434;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 414;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 433;
%                        Cntr = Cntr + 1;                                              
%                 end
%             elseif jjj == 2
%                 if XpsoP2(iii,is) == 0
%                     continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 423 322 431 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 424 422 331 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 423;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 424;
%                        Cntr = Cntr + 1;                                                                                            
%                 end
%             end
%         elseif iii == 5
%             if jjj == 1
%                 if XpsoP1(iii,is) == 0
%                     continue
%                 elseif XpsoP1(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 513 412 521 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 534 532 441 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 513;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 534;
%                        Cntr = Cntr + 1;                                                                     
%                 elseif XpsoP1(iii,is) == 2
%                        fprintf(fid,'element corotTrussSection 514 512 421 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 533 432 541 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 514;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 533;
%                        Cntr = Cntr + 1;                                                                                            
%                 else
%                        fprintf(fid,'element corotTrussSection 513 412 521 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 534 532 441 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 514 512 421 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 533 432 541 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 513;                       
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 534;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 514;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 533;
%                        Cntr = Cntr + 1;                                                                  
%                 end
%             elseif jjj == 2
%                 if XpsoP2(iii,is) == 0
%                     continue
%                 elseif XpsoP2(iii,is) == 1
%                        fprintf(fid,'element corotTrussSection 523 422 531 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        fprintf(fid,'element corotTrussSection 524 522 431 %d 4\r\n',    BrSection(XpsoBr(iii,is),1));
%                        CorE(Cntr,1) = 523;
%                        Cntr = Cntr + 1;
%                        CorE(Cntr,1) = 524;
%                        Cntr = Cntr + 1;                                                                                            
%                        
%                 end
%             end
%         end
%     end
% end

CorEE = zeros(Cntr-1,4,'gpuArray');

for ii = 1:Cntr-1            
    CorEE(ii,1) = CorE(ii,1);
    CorEE(ii,2) = CorE(ii,2);
    CorEE(ii,3) = CorE(ii,3);
    CorEE(ii,4) = CorE(ii,4);
end          
                   
WF = Gamma*(Ls*(4*(CSection(XpsoC(1,is),5) + CSection(XpsoC(2,is),5) + CSection(XpsoC(3,is),5) + CSection(XpsoC(4,is),5) + CSection(XpsoC(5,is),5) + ...
                   CSection(XpsoC(6,is),5) + CSection(XpsoC(7,is),5) + CSection(XpsoC(9,is),5) + CSection(XpsoC(10,is),5) + CSection(XpsoC(11,is),5) + ...
                   CSection(XpsoC(12,is),5) + CSection(XpsoC(13,is),5) + CSection(XpsoC(14,is),5) + CSection(XpsoC(15,is),5) + CSection(XpsoC(6,is),5) + ...
                   CSection(XpsoC(7,is),5) + CSection(XpsoC(9,is),5) + CSection(XpsoC(10,is),5)) + 2*(CSection(XpsoC(8,is),5) + CSection(XpsoC(16,is),5)))) ...
+ Gamma*(Lb*(3*NS*(BSection(XpsoB(1,is),1))));
                   
for iii = 1:size(CorEE,1)
    WF = WF + Gamma*Lt*CorEE(iii,4);
end   

                       


fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

% fprintf(fid,'set omega2list [eigen 1] \r\n');
% fprintf(fid,'set omega2 [lindex $omega2List 0] \r\n');
% fprintf(fid,'set Tperiod [expr 2*3.1415/sqrt($omega2)] \r\n');

fprintf (fid,  'pattern Plain 1 Linear { \r\n');

Load = gpuArray(-(1.2*2500+1.6*1000));

for ii = 1:size(BE)
    fprintf(fid,'eleLoad -ele %d -type -beamUniform %d \r\n', BE(ii,1), Load);
end

fprintf (fid,  '} \n\n');

%Define Recorder

% fprintf(fid,'recorder Element -file ElementResults.txt -time ele %d %d %d force \r\n', 1, size(CElements,1));
% fprintf(fid,'recorder Node -file Eigen1.txt node 11 12 13 14 15 -dof 1 "eigen 1" \r\n');

fprintf (fid,'recorder Element -file BMoment.out -ele 112 localForce \r\n');
fprintf(fid,'recorder Element -file CForce.out -ele');

Cntr = 1;

for iii = 1:nvc/2
    fprintf(fid,' %d',CE(Cntr,1));
    fprintf(fid,' %d',CE(Cntr + 1,1));
    Cntr = Cntr + 2*(NB+1);
end

fprintf(fid,' localForce \r\n');

% fprintf (fid,'recorder Element -file CForce.out -ele 111 311 511 121 321 521 localForce \r\n');

fprintf (fid,'recorder Element -file BrForce.out -ele');

Cntr = 1;

for iii = 1:size(CorEE,1)
    fprintf(fid,' %d',CorEE(Cntr,1));
    Cntr = Cntr + 1;
end
fprintf(fid,' axialForce \r\n');


%Gravity Analysis
fprintf (fid,'constraints Plain\r\n');
fprintf (fid,'numberer Plain\r\n');
fprintf (fid,'system SparseGEN \r\n');
fprintf (fid,'test NormDispIncr 1.e-6 10\r\n');
fprintf (fid,'algorithm Newton  \r\n');
fprintf (fid,'integrator LoadControl 1 \r\n');
fprintf (fid,'analysis Static \r\n');
fprintf (fid,'analyze 1 \r\n');
fprintf(fid,'exit \r\n');



fclose(fid);	