function [x_final,MSE] = Deblurring_GSR(i,path,y,Opts)

mu = Opts.mu;
x_org = Opts.org;
IterNums = Opts.IterNums;
initial = Opts.initial;

x = initial;

c = zeros(size(y));

MSE = zeros(1,IterNums+1);

MSE(1) = sum(sum((x-x_org).^2))/numel(x);

fprintf('Initial PSNR = %f\n',csnr(x,x_org,0,0));

result=fopen(strcat(path,'\data.txt'),'a+');
for Outloop = 1:IterNums
    
    w = GSR_Solver_Deblur(x-c,Opts);
    
    r = y + mu*(w+c);
    x = r/(mu+1);
    
    c = c + (w - x);
    
    x_resid = x - x_org;
    MSE(Outloop+1) =  (x_resid(:)'*x_resid(:))/numel(x);
    
    Final_Name_GSR = strcat('_IterNum_',num2str(Outloop),'_PSNR_',num2str(csnr(x,x_org,0,0)));
    fprintf('iter number = %d, PSNR = %f\n',Outloop,csnr(x,x_org,0,0));
    fprintf(result,'IterNum = %d, PSNR = %0.5f\r\n',Outloop,csnr(x,x_org,0,0));
end

x_final = x;

fprintf('Final PSNR = %f\n',csnr(x_final,x_org,0,0));

fprintf(result,'\r\n\r\n');
fclose('all');
if mod(i,10)==0
    save(strcat(path,'\',Final_Name_GSR,'.mat'),'x_final');
end
end



