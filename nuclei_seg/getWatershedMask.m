function M=getWatershedMask( I )
%GETWATERSHEDMASK Uses the fast veta watershed approach (Cheng Lu's) to segment
%nuclei of an image and saves the corresponding mask.

%addpath(genpath('nuclei_seg/veta_watershed'));
%addpath(genpath('nuclei_seg/GeneralLoG'));
%addpath(genpath('nuclei_seg/staining_normalization'));

%imgFile=[folder '/' image];
%curIM=imread(imgFile);

%[~, filename, ~] = fileparts(imgFile);

[w,h,~]=size(I);
[~,normI,~] = normalizeStaining(I,220,0.15);
%normRed=normI(:,:,1);
normRed=rgb2gray(normI);
%% using multi resolution watershed, speed up veta

p.scales=[4:2:14];
%disp('begin nuclei segmentation using watershed');
[nuclei, ~] = nucleiSegmentationV2(normRed,p);

M=zeros(w,h);

%imshow(curIM);
%hold on;
for k = 1:length(nuclei)
    nuc=nuclei{k};
    numPix=length(nuc);
    for ind=1:numPix
        M(nuc(ind,1),nuc(ind,2))=255;
    end
%    plot(nuclei{k}(:,2), nuclei{k}(:,1), 'g-', 'LineWidth', 2);
end
%hold off;

M = imfill(M,'holes');
%imwrite(mask,[folder '/' filename '_mask.png']);



end

