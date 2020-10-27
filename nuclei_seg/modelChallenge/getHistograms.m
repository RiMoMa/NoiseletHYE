function imgHist = getHistograms(img,msk)
% Get the local and contextual features of each Nuclei
[centroids,localFeatures,localLabels, contextFeatures, contextualLabels] = getNuclearFeatures(img,msk);
% Put together the both set of features
img_feat = [localFeatures contextFeatures];
% Convert any complex number, NaN and Inf to 0
img_feat(imag(img_feat) ~= 0) = 0;
img_feat(isnan(img_feat)) = 0;
img_feat(isinf(img_feat)) = 0;
% Size of the features
[~,feat_list] = size(img_feat);
% Size of the histogram's bins
bin_feat_clas = 10;
% Generate the histograms
imgHist = [];
for j =1:feat_list
    % Extract each feature
    img_tmp = img_feat(:,j);
    % Calculate the min and max value
    img_min_tmp = min(img_tmp);
    img_max_tmp = max(img_tmp);
    % Take the max and min and normalize
    img_stnd = (img_tmp - img_min_tmp)/(img_max_tmp - img_min_tmp);
    img_histo = histogram(img_stnd, bin_feat_clas, 'BinLimits',[0,1]);
    img_mom = img_histo.Values;
    imgHist = [imgHist img_mom];
end
end