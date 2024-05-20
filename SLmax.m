function SLmax_output = SLmax(xM1,yM2)
% Shifted Lmax 
% Returns the lexicographic maximum of 2 vectors xM1&yM2, left-shifted by
% one element
% Vectors xM1&yM2 must have the same size
 SLmax_output = xM1(2:end);
 for k=1:size(xM1,2)
    if xM1(k)>yM2(k)
       SLmax_output = xM1(2:end);
       break
    elseif xM1(k)<yM2(k)
       SLmax_output = yM2(2:end);
       break
    end
 end
end
