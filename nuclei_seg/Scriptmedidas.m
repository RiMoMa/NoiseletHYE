path_im = 'A12_00Dc.tiff';

im = imread(path_im);

M=getWatershedMask( im,'minScale',5,8);%mascar binaria

imshow(M)
title('mascara de nucleos')
 %%%sobre la mascara M se pueden realizar medidas
M= im2bw(M);

nuclei_centroids = regionprops(M,'centroid');
area_nuclei =regionprops(M,'Area');

%%%% para ver que nucleo es en la imagen cada medida

L = bwlabeln(M);
imagesc(L)