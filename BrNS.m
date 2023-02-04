nBr1 = 0;
nBr2 = 0;
nBr3 = 0;
nBr4 = 0;
nBr5 = 0;

for iii = 1:size(CorEE,1)
    if 0 < CorEE(iii,1) - 100 < 100
       nBr1 = nBr1 + 1;
    elseif 100 < CorEE(iii,1) - 100 < 200
       nBr2 = nBr2 + 1;
    elseif 200 < CorEE(iii,1) - 100 < 300
       nBr3 = nBr3 + 1;
    elseif 300 < CorEE(iii,1) - 100 < 400
       nBr4 = nBr4 + 1;
    elseif 400 < CorEE(iii,1) - 100 < 500
       nBr5 = nBr5 + 1;
    end
end

