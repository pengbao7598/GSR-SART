function [ trunVec ] = trunFuncTh( inVec, lowThd, upThd )
%TRUNFUNCTH truncate vector with lower and upper threshold 
%   inVec    input vector
%   lowThd   lower threshold
%   upThd    upper threshold
%   trunVec  truncated vector
temp = inVec;
temp(find(temp<=lowThd)) = [];
temp(find(temp>=upThd)) = [];
trunVec = temp;

end

