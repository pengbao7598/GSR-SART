function [ systemMatrix ] = generateParalSysMatrix( h1,h2,d,num,length,width,M,N,angNum  )
%GENERATEPARALSYSMATRIX Generate system matrix arrording to the geometrical
%                     configuration
%   type parallelBeam or fanBeam
%   h1 distance from the source to the center of rotation
%   h2 distance from the center of rotation to the detector panel
%   d  length of the detector panel
%   num   number of bins on the detector
%   length  length of image (cm)
%   width   width of image  (cm)
%   [M N] is size of image (pixels)
%   angNum sampling number of angles

scale1 = width/M;
scale2 = length/N;
if scale1 == scale2
    scale = scale1;
else
    error('The scale is not coherent! Please check!');
end

yPlane = -width/2:width/M:width/2;
xPlane = -length/2:length/N:length/2;
projNum = angNum*num;
systemMatrix = sparse(projNum,M*N);
rowPointer = 1;
% flag=1;
for angTemp = 0:pi/angNum:pi-pi/angNum
    [sourPos,endPos] = generateParallelIndices(h1,h2,d,num,angTemp);
    
    for endNum = 1:num
        sourceSinglePos=sourPos(endNum,:);
        endSinglePos = endPos(endNum,:);
        [~,~,~,systemMatrixRow] = extractRadiPathSiddon(sourceSinglePos,endSinglePos,xPlane,yPlane);
        systemMatrix(rowPointer,:) = systemMatrixRow;
        rowPointer = rowPointer+1;
    end
%     flag=flag+1
end

end

