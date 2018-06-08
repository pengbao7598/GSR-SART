function [ filteredVec ] = filterVector( inVec )
%FILTERVECTOR filter the input vector -- remove the items approximately equal
%   inVec   input vector
%   filteredVec   filtered vector
[m n] = size(inVec);
temp = [];
if (m>1) && (n>1)
    error('the dimension error in function filterVector()');
end
if (m>1) && (n==1)
    inVec = inVec';
end
for i = 1:(n-1)
    if abs(inVec(1,i+1)-inVec(1,i))<10^(-7)
        temp(end+1) = inVec(i);
    end
end
numRep = size(temp,2);
for j = 1:numRep
    inVec(find(inVec==temp(j))) = [];
end
filteredVec = inVec;
end

