function [ x ] = sart( A, w, y, x0, t)
%SART Summary of this function goes here
%   Detailed explanation goes here
A_i = sum(A,2);  % sum of each row
A_j = sum(A,1)';  % sum of each column
y_temp = A*x0;
AT=A';

if t == 1
    x = x0+w.*(sum(AT*((y-y_temp)./A_i),1)./A_j)';
else
    x_temp = x0;
    for z = 1:t
        %         a=(y-y_temp);
        %         b=a./A_i;
        %         c=AT*b;
        %         d=c/A_j;
        %         x=x_temp+w*d;
        x = x_temp+w*(AT*((y-y_temp)./A_i)./A_j);
        y_temp = A*x;
        x_temp = x;
        x_temp = max(x_temp,0);
    end
end
x= reshape(x,256,256);
end



