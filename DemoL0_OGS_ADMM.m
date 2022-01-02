
% 7/9/2019
% This is the demo file to run image denoising and deblurring with impulse
% noise. The algorithm is the L0-OGSTV.

% L0 --> Data fidelity term from Yuan and Ghanem (for impulse noise)
% OGSTV --> Regularization term (supress staircase artifacts)

% See the function "L0_OGS_ADMM.m" for further details on L1-OGSTV


%function DemoL0_OGS_ADMM(lam)


%{
clc;
clear all;
close all;
%}



NN = 10;
sumPsnr = 0;
sumSSIM =0;
sumSnr  = 0;
sumIt   = 0; % to calculate average itration to converge
sumTime = 0;

saveRes = 'no';

%for j = 1:NN

imageName = 'boats.bmp';    

Img = imread(imageName);

if size(Img,3) > 1
    Img = rgb2gray(Img);
end

[row, col] = size(Img);

row = int2str(row);
col = int2str(col);

imageSize = [row 'x' col];

K = fspecial('gaussian', [7 7], 5); % Gaussian Blur
%K     =   fspecial('average',7); % For denoising
f1 = imfilter(Img,K,'circular');
f1 = double(f1);

%f = imnoise(f,'salt & pepper',0.5);

f  = impulsenoise(f1,0.5,0);
f = double(f);

O = ones(size(Img));
O(f == 255) = 0;
O(f == 0) = 0;

Img = double(Img)/255;
f = f/255;


opts.lam       = 0.11; % for denoising try starting with these  
opts.grpSz     = 3; % OGS group size
opts.Nit       = 1000;
opts.Nit_inner = 5;
opts.tol       = 1e-4;
opts.O        = O;

% main function

    
    out = L0_OGS_ADMM(f, Img, K, opts);
   
   %{
    sumPsnr = sumPsnr + out.psnrRes;
    sumSSIM = sumSSIM + out.ssimRes;
    sumIt   = sumIt   + out.OverallItration ;
    sumSnr  = sumSnr  + out.snrRes;
    sumTime = sumTime + out.cpuTime;
    
    
    clear Img out f opts N K;
    
    
    
end
    avePsnr = sumPsnr/NN;
    aveSSIM = sumSSIM/NN;
    aveIt   = sumIt/NN; 
    aveTime = sumTime/NN;
    aveSnr  = sumSnr/NN;
%%
%}
%
figure;
imshow(out.sol,[]),
title(sprintf('ogs2d\\_tv (PSNR = %3.2f dB,SSIM = %3.3f, cputime %.2f s) ',...
                      psnr_fun(out.sol*255,Img*255),ssim_index(out.sol*255,Img*255)))
 
 %
figure;                   
imshow(Img,[])
title('Original Image');

figure;
imshow(f,[])
title('Blur + Noisy');
%}

%}                
                   
                   %Some result options

fileName = 'boats_result_0_9snp_denoising_15x15gaussion_5_groupsize_3_L0_OGS_admm.txt';
        
    switch saveRes
    case 'yes'
        %mkdir('Result');
        
        if exist(fileName,'file')
                  fileID = fopen(fileName,'a');
                  fprintf(fileID,'%6s   %9s     %9.4f      %9.4f      %9.4f      %9.2f      %9.4f \n',imageName,imageSize, psnr_fun(out.sol*255,Img*255),ssim_index(out.sol*255,Img*255),out.cpuTime,out.OverallItration, opts.lam );
                  fclose(fileID);
        else
        
                fileID = fopen(fileName,'w');
                fprintf(fileID,sprintf('%6s     %9s    %9s    %9s     %9s        %9s        %9s  \n\n','Image','Size','PSNR','SSIM','Time','Iteration','lam'));
                fprintf(fileID,'\n\n');
                fprintf(fileID,'%6s   %9s    %9.4f      %9.4f    %9.4f   %9.2f    %9.4f \n',imageName, imageSize, psnr_fun(out.sol*255,Img*255),ssim_index(out.sol*255,Img*255),out.cpuTime,out.OverallItration, opts.lam);
                fclose(fileID);
        end
        
    case 'no'
        %break;
         return;
        
    otherwise
        warning('Incorrect options. Options are "yes" or "no"');
       % break;
        return;
       
end

