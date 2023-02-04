
if iter == 1
    Gb = min(SOFW);
    locg = find(SOFW == Gb);
    for iii = 1:np
        XgbestC(:,iii) = XpsoC(:,locg(1));
        XgbestB(:,iii) = XpsoB(:,locg(1));
        XgbestBr(:,iii) = XpsoBr(:,locg(1));
        XgbestP1(:,iii) = XpsoP1(:,locg(1));
        XgbestP2(:,iii) = XpsoP2(:,locg(1));
    end
else
    minsupf = min(SOFW);
    locg2 = find(SOFW == minsupf);
    if minsupf < Gb
       Gb = minsupf;
       for iii = 1:np
           XgbestC(:,iii) = XpsoC(:,locg2(1));
           XgbestB(:,iii) = XpsoB(:,locg2(1));
           XgbestBr(:,iii) = XpsoBr(:,locg2(1));
           XgbestP1(:,iii) = XpsoP1(:,locg2(1));
           XgbestP2(:,iii) = XpsoP2(:,locg2(1));
       end
    end
end


% if iter == 1
%     Gb = min(SOFW);
%     locg = find(SOFW == Gb);
%     for iii = 1:np
%         XgbestB(:,iii) = XpsoB(:,locg(1));
%     end
% else
%     minsupf = min(SOFW);
%     locg2 = find(SOFW == minsupf);
%     if minsupf < Gb
%        Gb = minsupf;
%        for iii = 1:np
%        end
%     end
% end

% 
% if iter == 1
%     Gb = min(SOFW);
%     locg = find(SOFW == Gb);
%     for iii = 1:np
%         XgbestBr(:,iii) = XpsoBr(:,locg(1));
%     end
% else
%     minsupf = min(SOFW);
%     locg2 = find(SOFW == minsupf);
%     if minsupf < Gb
%        Gb = minsupf;
%        for iii = 1:np
%            XgbestBr(:,iii) = XpsoBr(:,locg2(1));
%        end
%     end
% end
% 
% 
% if iter == 1
%     Gb = min(SOFW);
%     locg = find(SOFW == Gb);
%     for iii = 1:np
%         XgbestP1(:,iii) = XpsoP1(:,locg(1));
%     end
% else
%     minsupf = min(SOFW);
%     locg2 = find(SOFW == minsupf);
%     if minsupf < Gb
%        Gb = minsupf;
%        for iii = 1:np
%            XgbestP1(:,iii) = XpsoP1(:,locg2(1));
%        end
%     end
% end
% 
% 
% if iter == 1
%     Gb = min(SOFW);
%     locg = find(SOFW == Gb);
%     for iii = 1:np
%         XgbestP2(:,iii) = XpsoP2(:,locg(1));
%     end
% else
%     minsupf = min(SOFW);
%     locg2 = find(SOFW == minsupf);
%     if minsupf < Gb
%        Gb = minsupf;
%        for iii = 1:np
%            XgbestP2(:,iii) = XpsoP2(:,locg2(1));
%        end
%     end
% end
% 
Xgbest = [XgbestC;XgbestB;XgbestBr;XgbestP1;XgbestP2];