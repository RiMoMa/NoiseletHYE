% DEMO for extracting nuclei local and contextual information from an
% image and generate histograms of the extracted features
% 2018
imgDir = 'image/';
maskDir = 'mask/';
% Load the image
img = imread([imgDir 'n001.png']); 
% Load the mask
%msk = imread([maskDir 'n001_mask.png']);
% Create the mask (Watershed-Veta)
msk = getWatershedMask(img,true,4,10);
% Generate the histograms of the image
imgH = getHistograms(img,msk);
% Generate the celularity Features
imgC = getCellularityImg(msk);
% Generate the Global Zernikes Features
%imgZ = add zernike global code
% Put the extracted Features together
img_features = [imgH imgC]; %imgZ
% Normalize the dataset
img_features_norm = normc(img_features); %imgZ