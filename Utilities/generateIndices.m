function [ sourPos,endPos ] = generateIndices( h1,h2,d,num,theta )
%GENERATEINDICES Generate the indice couples of fan beam projection
%   h1 distance from the source to the center of rotation
%   h2 distance from the center of rotation to the detector panel
%   d  length of the detector panel
%   num   number of bins on the detector
%   theta  rotation angle of source
s_xPos = -cos(theta)*h1;
s_yPos = sin(theta)*h1;
if abs(s_xPos)<10^(-7)
    s_xPos = 0;
end

if abs(s_yPos)<10^(-7)
    s_yPos = 0;
end

sourPos = [s_xPos s_yPos];        % calculate the coordinate of source
endPos = zeros(num,2);      % store the coordinates of end points into a two-dimensional matrix
d_interval = d/(num-1);          % interval of the detecors
ite = 1;

for t = d/2:-d_interval:-d/2
    ang = atan(t/h2);
    dis = (t^2+h2^2)^0.5;
    xPos = dis*cos(ang-theta);
    yPos = dis*sin(ang-theta);
    if abs(xPos)<10^(-7)
        xPos = 0;
    end
    if abs(yPos)<10^(-7)
        yPos = 0;
    end
    endPos(ite,:) = [xPos yPos];
    ite = ite+1;
end

end

