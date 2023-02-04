
if iter == 1
    fitness = SOFW;
    XpbestC = XpsoC;
    XpbestB = XpsoB;
    XpbestBr = XpsoBr;
    XpbestP1 = XpsoP1;
    XpbestP2 = XpsoP2;
else
    for ipb = 1:np
        if SOFW(ipb) <= fitness(ipb)
            fitness(ipb) = SOFW(ipb);
            XpbestC(:,ipb) = XpsoC(:,ipb);
            XpbestB(:,ipb) = XpsoB(:,ipb);
            XpbestBr(:,ipb) = XpsoBr(:,ipb);
            XpbestP1(:,ipb) = XpsoP1(:,ipb);
            XpbestP2(:,ipb) = XpsoP2(:,ipb);
        end
    end
end

% if iter == 1
%     fitness = SOFW;
%     XpbestB = XpsoB;
% else
%     for ipb = 1:np
%         if SOFW <= fitness(ipb)
%             fitness(ipb) = SOFW(ipb);
%             XpbestB(:,ipb) = XpsoB(:,ipb);
%         end
%     end
% end

% if iter == 1
%     fitness = SOFW;
%     XpbestBr = XpsoBr;
% else
%     for ipb = 1:np
%         if SOFW <= fitness(ipb)
%             fitness(ipb) = SOFW(ipb);
%             XpbestBr(:,ipb) = XpsoBr(:,ipb);
%         end
%     end
% end

% if iter == 1
%     fitness = SOFW;
%     XpbestP1 = XpsoP1;
% else
%     for ipb = 1:np
%         if SOFW <= fitness(ipb)
%             fitness(ipb) = SOFW(ipb);
%             XpbestP1(:,ipb) = XpsoP1(:,ipb);
%         end
%     end
% end

% if iter == 1
%     fitness = SOFW;
%     XpbestP2 = XpsoP2;
% else
%     for ipb = 1:np
%         if SOFW <= fitness(ipb)
%             fitness(ipb) = SOFW(ipb);
%             XpbestP2(:,ipb) = XpsoP2(:,ipb);
%         end
%     end
% end
