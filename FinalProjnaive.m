clear;
clc;
close all;

%%% Declaring parameters for the retargeting
minImgSize = 30;                % lowest scale resolution size for min(w, h)
outSizeFactor = [1, 0.65];		% the ration of output image
numScales = 5;                 % number of scales (distributed logarithmically)

%% Preparing data for the retargeting
image = imread('SimakovFarmer.png');
[h, w, ~] = size(image);

targetSize = round(outSizeFactor .* [h, w]);
imageLab = uint8(image);%rgb2lab(image); % Convert the source and target Images
%imageLab = double(imageLab);%/255;

% Gradual Scaling - iteratively icrease the relative resizing scale between the input and
% the output (working in the coarse level).
%% STEP 1 - do the retargeting at the coarsest level
tic;
reducedRatio = [35, 47];
numScales1 = 10;
rcols = 10^((1/numScales1)*log10(minImgSize/reducedRatio(2)));

A = imresize(imageLab, reducedRatio);
B = A;
num =1;
for num = 1:numScales1
    T = A;
    B = imresize(B,[35,round(47*rcols^num)]);
    [B, ann, bnn] = search_vote_funcnaive(T,B,4);
end
disp(['Patch-Match time: ', num2str(toc), ' sec']);

 subplot(1,3,2); imshow(A);
 subplot(1,3,3); imshow(B);
 subplot(1,3,1); imshow(imageLab);
%{
%% STEP 2 - do resolution refinment 
% (upsample for numScales times to get the desired resolution)
redRatio = size(B);
rrows2 = 10^((1/numScales)*log10( targetSize(1)/redRatio(1) ));
rcols2 = 10^((1/numScales)*log10( targetSize(2)/redRatio(2) ));

rrowsact = 10^((1/numScales)*log10( h/redRatio(1) ));
rcolsact = 10^((1/numScales)*log10( w/redRatio(2) ));
for num = 1:numScales
    m = redRatio(1)*rrows2^num; n = redRatio(2)*rcols2^num;
    m2 = m; n2 = redRatio(2)*rcolsact^num;
    T = imresize(imageLab, [m2,n2]); 
    B = imresize(B, [m,n]);
    B = search_vote_func(T, B, 2); 
end


%% STEP 3 - do final scale iterations
% (refine the result at the last scale)
F = search_vote_func(imageLab,B, 3);
disp(['Patch-Match time: ', num2str(toc), ' sec']);

figure
imshow(F);
% %% Seam Carving
% tic;
% img = imageLab;
% Img_d = double(img);%/255;
% EMap = myEnergyFunc(Img_d);
% % figure,
% % imshow(EMap)
% % title('Energy Map')
% rC = w - targetSize(2);
% rR = h - targetSize(1);
% 
% rImg = mySeamCarveResize(Img_d,rC,rR);
% disp(['Seam Carving time: ', num2str(toc), ' sec']);
% 
% %%
% %imshow(rImg)
% C = [imageLab,rImg]*255;
% C = [F,C];
% figure
% imshow(C)
%}