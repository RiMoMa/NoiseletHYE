% DEMO 2 for extracting nuclei local and contextual information from all
% the test images of the Challange and generate histograms of the 
% extracted features
% 2018

% Paths
imgsDir = 'png_test/';
% Only if the mass are available
maskDir = 'png_test_masks/';
saveDir = 'png_test_feats/';

% Label for the type or level of cancer used
% 1 b = benign
% 2 is = insitu 
% 3 iv = invasive 
% 4 n = normal

% Load the images
gathimgs = dir([imgsDir '*png']);
for a = 1:length(gathimgs)
    imgp = [imgsDir gathimgs(a).name];
    [~, name, ~] = fileparts(imgp);
    mskp = [maskDir name '_mask.png'];
    disp(['Extracting Features from Image ' name]);
    % Load image and mask
    imgl = imread(imgp);
    mskl = imread(mskp);
    % If masks are not pressent, use watershed-veta. Save mask if need it
    %mskl = getWatershedMask(imgl,true,4,10);
    %imwrite(maskl,[maskDir name '_mask.png']);
    % get Histograms
    imgHistograms = getHistograms(imgl,mskl);
    % get Celularity
    imgCelularity = getCellularityImg(mskl);
    
    % get Global Zernikes
    [zer_features] = zernike_img_extraction(im2single(imgl)); 

    % imgGlobalZernike = add global zernikes code
    save([saveDir name '.mat'],'imgHistograms','imgCelularity','zer_features'); % 'imgGlobalZernike'
end