%   These MATLAB programs implement the Few-view CT reconstruction method as described in paper:
%   
%     Title: Few©\view CT reconstruction with group©\sparsity regularization
%     Author: Peng Bao, Jiliu Zhou, Yi Zhang
%   
% -------------------------------------------------------------------------------------------------------
% The software implemented by MatLab 2017a are included in this package.

% ------------------------------------------------------------------
% Author: Peng Bao
% Email:  pengbao7598@gmail.com
% Last modified by Peng Bao, June 2018

clear
clc
cur = cd;
addpath(genpath(cur));


ImgNo=1;

switch ImgNo
    case 1
        imageName = 'abdominal_image';
    case 2
        imageName = 'pelvic_image';
    case 3
        imageName = 'thoracic_image';
end

load(strcat('Test_Images/',imageName,'.mat'));

% figure,imshow(original_image,[850/3000 1250/3000]);

[row,col] = size(original_image);

Opts = [];
Opts.nangle=64;
Opts.nbin=512;

%generate system matrix
% SystemMatrix=Parallel_BeamSystemMatrix(30,150,col,row,1,1);
% SystemMatrix=generateSystemMatrix(80,80,150,70,row,col,row,col,30);
% SystemMatrix=generateSystemMatrix(100,100,300,Opts.nbin,row,col,row,col,Opts.nangle);   %96*96
% SystemMatrix=generateSystemMatrix(100,100,400,Opts.nbin,row,col,row,col,Opts.nangle);
% SystemMatrix=generateParalSysMatrix(200,200,200,Opts.nbin,row,col,row,col,Opts.nangle);
% load('128_60_150_SystemMatrix.mat')

proj = 1;
resolution = 256;
switch proj
    case 1
        Projection_method = 'Fan-Beam';
%         load('128_64_150_FanSysMat.mat')
%         load('128_64_150_FanSysMat2.mat')
%         load('128_64_256_FanSysMat2.mat');
        load('Data/256_64_512_FanSysMat2.mat');
%         load('256_30_512_FanSysMat2.mat');
        %         load('128_30_250_FanSysMat.mat')
    case 2
        Projection_method = 'Paral-Beam';
        % load('128_30_150_ParalSysMat.mat')
end

Opts.Phi = SystemMatrix;
Opts.row = row;
Opts.col = col;

% Get Measurements
y = SystemMatrix * original_image(:);

if ~isfield(Opts,'IterNum')
    %Opts.IterNum = 120;
    Opts.IterNums = 3;
end

if ~isfield(Opts,'mu')
    Opts.mu = 0.1;
end

if ~isfield(Opts,'lambda')
    Opts.lambda = 0.00005;
end

if ~isfield(Opts,'Inloop')
    Opts.Inloop = 50;
end

if ~isfield(Opts,'PatchSize')
    Opts.PatchSize = 8;
end

if ~isfield(Opts,'ArrayNo')
    Opts.ArrayNo = 40;
end

if ~isfield(Opts,'SlidingDis')
    Opts.SlidingDis = 4;
    Opts.Factor = 240*Opts.ArrayNo/60;
end

if ~isfield(Opts,'SearchWin')
    Opts.SearchWin = 20;
end

%params for sart
w=1;
t=500;

%main
% for j = 1 : 4
%     Opts.lambda = 0.05 * 10.^(j-1);
%     for k = 1 : 3
%         Opts.mu = 0.005 * 10.^(k-1);
%         Opts.SearchWin = 50+10*(k-1);
        
        Document_Name=strcat(Projection_method,'_',num2str(resolution),'_',imageName,'_lambda_',num2str(Opts.lambda),'_mu_',num2str(Opts.mu),'_nangle_',num2str(Opts.nangle),'_bin_',num2str(Opts.nbin),'_ArrayNo_',num2str(Opts.ArrayNo),'_SearchWin_',num2str(Opts.SearchWin));
        path=strcat('GSR_Result\',Document_Name);
        if ~exist(path)
            mkdir(path);
        end
        
        if exist(strcat(path,'\data.txt'),'file')
           delete( strcat(path,'\data.txt'));
        end
        
        for i=1:300
            result=fopen(strcat(path,'\data.txt'),'a+');      
            fprintf(result,'loop = %d \r\n',i);
            fprintf(result,'lambda = %d, mu = %0.5f\r\n',Opts.lambda,Opts.mu);
            fprintf('Loop:%d\n',i);
            %SART
            if(i==1)
                x0=zeros(row*col,1);
            else
                x0=reconstructed_image(:);
            end
            reconstructed_image_sart = sart( SystemMatrix, w, y, x0, t);
            
%             fprintf('RMSE = %0.5f',RMS(original_image,reconstructed_image_sart));
%             fprintf('RMSE = %0.5f',psnr(original_image,reconstructed_image_sart));
%             figure;imshow(reconstructed_image_sart,[0 0.05]);

            if mod(i,10)==0
                save(strcat(path,'\sart_PSNR_',num2str(csnr(original_image,reconstructed_image_sart,0,0)),'.mat'),'reconstructed_image_sart');
            end
            
            Opts.initial = double(reconstructed_image_sart);
            
            Opts.org = original_image;
            
            fprintf('Initial PSNR = %0.2f\n',csnr(Opts.org,Opts.initial,0,0));
            
            % Invoke Proposed GSR Alogorithm for Regularization
            disp('Beginning of GSR Algorithm for Regularization');
            
            reconstructed_image=Deblurring_GSR(i,path,reconstructed_image_sart,Opts);
            psnr = PSNR(original_image, reconstructed_image);
            %figure(2); imshow(double(reconstructed_image)); title('reconstructed_image');
        end
        
        fclose('all');
%     end
% end
disp('End of GSR');

