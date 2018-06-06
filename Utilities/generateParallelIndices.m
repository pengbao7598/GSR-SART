function [ sourPos,endPos ] = generateParallelIndices( h1,h2,d,num,theta )
%GENERATEINDICES Generate the indice couples of fan beam projection
%   h1 dis1tance from the source to the center of rotation
%   h2 dis1tance from the center of rotation to the detector panel
%   d  length of the detector panel
%   num   number of bins on the detector
%   theta  rotation angle of source

sourPos = zeros(num,2);        % calculate the coordinate of source
endPos = zeros(num,2);      % store the coordinates of end points into a two-dimensional matrix
d_interval = d/(num-1);          % interval of the detecors
ite = 1;

for t = d/2:-d_interval:-d/2
    ang1=atan(t/h1);
    dis1=(t^2+h1^2)^0.5;
    sour_xPos=-dis1*cos(ang1+theta);
    sour_yPos=dis1*sin(ang1+theta);
    
    if abs(sour_xPos)<10^(-7)
        sour_xPos = 0;
    end
    if abs(sour_yPos)<10^(-7)
        sour_yPos = 0;
    end
    sourPos(ite,:) = [sour_xPos sour_yPos];
    
    ang2 = atan(t/h2);
    dis2 = (t^2+h2^2)^0.5;
    end_xPos = dis2*cos(ang2-theta);
    end_yPos = dis2*sin(ang2-theta);
   
    if abs(end_xPos)<10^(-7)
        end_xPos = 0;
    end
    if abs(end_yPos)<10^(-7)
        end_yPos = 0;
    end
    endPos(ite,:) = [end_xPos end_yPos];
    ite = ite+1;
end

end


