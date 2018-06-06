function [ l,i,j,systemMatrixRow ] = extractRadiPathSiddon( posStart, posEnd, xPlane, yPlane )
%EXTRACTRADIPATHSIDDON Extract radiological path using Siddon algorithm
%   posStart  the indices of start point (x1,y1);
%   posEnd  the indices of end point (x2,y2);
%   xPlane  the parallel x planes of target object
%   yPlane  the parallel y planes of target object
%   l  length vector
%   i  x indices vector
%   j  y indices vector
%   systemMatrixRow one row of system matrix

[dim1,nx] = size(xPlane);           % nx is the length of xPlane
[dim2,ny] = size(yPlane);           % ny is the length of yPlane
x1 = posStart(1,1);
y1 = posStart(1,2);
x2 = posEnd(1,1);
y2 = posEnd(1,2);
dx = xPlane(1,2)-xPlane(1,1);
dy = yPlane(1,2)-yPlane(1,1);

%% normal situation
if (x1-x2 ~= 0)&&(y1-y2 ~= 0)
    alpha_x1 = (xPlane(1,1)-x1)/(x2-x1);
    alpha_xN = (xPlane(1,nx)-x1)/(x2-x1);         % refer to eqn.(4) in the reference
    alpha_y1 = (yPlane(1,1)-y1)/(y2-y1);
    alpha_yN = (yPlane(1,ny)-y1)/(y2-y1);
    
    systemMatrixRow = zeros(1,(ny-1)*(nx-1));
    i = 0;
    j = 0;
    l = 0;
    
    alpha_min = max([0 min(alpha_x1,alpha_xN) min(alpha_y1,alpha_yN)]);
    alpha_max = min([1 max(alpha_x1,alpha_xN) max(alpha_y1,alpha_yN)]);       % refer to eqn.(5) in the reference
    
    if alpha_min >= alpha_max
        return;
    end
    
    if (x2-x1)>=0
        i_min = nx-(xPlane(1,nx)-alpha_min*(x2-x1)-x1)/dx;
        i_max = 1+(x1+alpha_max*(x2-x1)-xPlane(1,1))/dx;
    else
        i_min = nx-(xPlane(1,nx)-alpha_max*(x2-x1)-x1)/dx;
        i_max = 1+(x1+alpha_min*(x2-x1)-xPlane(1,1))/dx;                    % refer to eqn.(6) in the reference x part
    end
    
    if (y2-y1)>=0
        j_min = ny-(yPlane(1,ny)-alpha_min*(y2-y1)-y1)/dy;
        j_max = 1+(y1+alpha_max*(y2-y1)-yPlane(1,1))/dy;
    else
        j_min = ny-(yPlane(1,ny)-alpha_max*(y2-y1)-y1)/dy;
        j_max = 1+(y1+alpha_min*(y2-y1)-yPlane(1,1))/dy;                    % refer to eqn.(6) in the reference y part
    end
    
    i_min = ceil(i_min);                % special part for some certain situation  2012-12-18 0:11
    i_max = floor(i_max);
    
    j_min = ceil(j_min);
    j_max = floor(j_max);
    
    alpha_x = [];
    alpha_y = [];
    
    if (x2-x1)>0
        for t = 1:(i_max-i_min+1)
            alpha_x(1,t) = (xPlane(1,i_min+t-1)-x1)/(x2-x1);               % refer to eqn.(7) in the reference x part
        end
    else
        for t = 1:(i_max-i_min+1)
            alpha_x(1,t) = (xPlane(1,i_max-t+1)-x1)/(x2-x1);
        end
    end
    
    if (y2-y1)>0
        for t = 1:(j_max-j_min+1)
            alpha_y(1,t) = (yPlane(1,j_min+t-1)-y1)/(y2-y1);               % refer to eqn.(7) in the reference y part
        end
    else
        for t = 1:(j_max-j_min+1)
            alpha_y(1,t) = (yPlane(1,j_max-t+1)-y1)/(y2-y1);
        end
    end
    alpha_x = trunFuncTh(alpha_x,alpha_min,alpha_max);
    alpha_y = trunFuncTh(alpha_y,alpha_min,alpha_max);
    alpha_temp = sort([alpha_x alpha_y],2);
    alpha = [alpha_min alpha_temp alpha_max];                              % refer to eqn.(8) in the reference
    
    alpha = unique(alpha);
    alpha = filterVector(alpha);
    [m n] = size(alpha);
    
    %    n = (i_max-i_min+1)+(j_max-j_min+1)+2;                                 % refer to eqn.(9) in the reference
    
    d_12 = ((x2-x1)^2+(y2-y1)^2)^0.5;                                      % refer to eqn.(11) in the reference
    
    for t = 1:(n-1)
        l(1,t) = (alpha(1,t+1)-alpha(1,t))*d_12;                           % refer to eqn.(10) in the reference
        if abs(l(1,t))<10^(-7)
            l(1,t) = 0;
        end
%         t
        alpha_mid = (alpha(1,t+1)+alpha(1,t))/2;                           % refer to eqn.(13) in the reference
        i(1,t) = 1+(x1+alpha_mid*(x2-x1)-xPlane(1,1))/dx;                  % refer to eqn.(12) in the reference x part
        j(1,t) = 1+(y1+alpha_mid*(y2-y1)-yPlane(1,1))/dy;                  % refer to eqn.(12) in the reference y part
        numPixel = floor(ny-j(1,t))*(nx-1)+floor(i(1,t));                  %calculate the numeber conrresponding to system matrix
        systemMatrixRow(1,numPixel) = l(1,t);                              %set value of system matrix
    end
end

%% when the radiological path is parallel with y axis
if (x2-x1 == 0)&&(y2-y1 ~= 0)
    systemMatrixRow = zeros(1,(ny-1)*(nx-1));
    if (x1<xPlane(1,1))||(x1>xPlane(1,nx))
        l = zeros(1,nx-1);    % length vector
        i = zeros(1,nx-1);    % x indices vector
        j = zeros(1,nx-1);    % y indices vector
    else
        i = x1;
        j = (yPlane(1,1)+dy/2):1:yPlane(1,ny)-dy/2;   %% problem here Yan Liu
        l = zeros(1,ny-1);
        for t = 1:(ny-1)
            numPixel = floor(yPlane(1,ny)-j(1,t))*(nx-1)+ceil(i-xPlane(1,1)); %% problem here Yan Liu
            systemMatrixRow(1,numPixel) = dy;
            l(1,t) = dy;
        end
    end
end

%% when the radiological path is parallel with x axis
if (x2-x1 ~= 0)&&(y2-y1 == 0)
    systemMatrixRow = zeros(1,(ny-1)*(nx-1));
    if (y1<yPlane(1,1))||(y1>yPlane(1,ny))
        l = zeros(1,ny-1);    % length vector
        i = zeros(1,ny-1);    % x indices vector
        j = zeros(1,ny-1);    % y indices vector
    else
        j = y1;
        i = (xPlane(1,1)+dx/2):1:xPlane(1,nx)-dx/2;   %% problem here Yan Liu
        l = zeros(1,nx-1);
        for t = 1:(nx-1)
            numPixel = floor(yPlane(1,ny)-j)*(nx-1)+ceil(i(1,t)-xPlane(1,1)); %% problem here Yan Liu
            systemMatrixRow(1,numPixel) = dx;
            l(1,t) = dx;
        end
    end
end

end

