function [labelsAll,histograms] = NoiseletsPLSAHistogramImg(im,Scales,WinPlsa,GroundT,vocab)
%%% Sacar Mascaras
addpath('PLSA_TEM_EM')
addpath('FNT')
%%% OPTIONS %%%
%K_clusters = 8;

%[Inorm H E] = normalizeStaining(im);
%[H,E] = ColorSepCimalab(im);
 %labim = rgb2lab(H);
 %ime = labim(:,:,2);
 %%ime = ((ime-min(ime(:)))/(max(ime(:))-min(ime(:))))*255;
 %ime=uint8(ime);
%ime = im(:,:,1);
% %%%%%%test 
% I = ime
% thresh = multithresh(I,2);
% valuesMax = [thresh max(I(:))];
% [quant8_I_max, index] = imquantize(I,thresh,valuesMax);
% valuesMin = [min(I(:)) thresh];
% quant8_I_min = valuesMin(index);
% imshowpair(quant8_I_min,quant8_I_max,'montage')
% title('Minimum Interval Value           Maximum Interval Value')
% ime = quant8_I_min;

win_size=WinPlsa;



[aa, bb, cc] = size(im);
andregions = ones(aa,bb);

    [Inorm,H,E] = normalizeStaining(im,im,220,0.06);
    %[H,E] = ColorSepCimalab(im);
    %labim = rgb2lab(H);
    ime = H(:,:,1);

    %%% Extract Noiselets from tiles
    [pila_patches_orig,coords,lienzo,pila_patches_gray,labels] = ExtractTilesAndNoiseletsLabels(Scales,ime,aa,bb,GroundT);
    
      %  dataset Monuseg
   % [XFeatures,imLabels] = NoiseletsFeaturesLabels(ime,Scales,WinPlsa,groundT);
    abs_pila = real(pila_patches_orig); %Se saca la Magnitud
    angle_pila = imag(pila_patches_orig); %SE SACA EL ANGULO, se le su

    %vector a clusterizacion
    X=zeros(size(pila_patches_orig,1),size(pila_patches_orig,2)*2);
    X(:,1:2:end-1)=abs_pila;
    X(:,2:2:end) = angle_pila; 
    XFeatures = X;

    %%histogram building
    [~,idx] = pdist2(vocab,XFeatures,'cityblock','Smallest',1);
    
    
  
       
       
       lienzo2 =lienzo


    
    %% Asigna la ocurrencia al tile
    for L = 1:length(coords)

        aux_point = coords(L,:);
        lienzo(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
            idx(L)*ones(win_size,win_size);
        
        
        lienzo2(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
            labels(L)*ones(win_size,win_size);

    end   
   
    
    
    %%%% PLSA Analysis
    
    
    [pila_patches_PLSA,coordsPLSA,lienzoPLSA,histograms] = ExtractTilesOneDPyramidal(WinPlsa,lienzo,aa,bb);
    [pila_patches_PLSA2,coordsPLSA2,lienzoPLSA2,histograms2] = ExtractTilesOneDPyramidal(WinPlsa,lienzo2,aa,bb);
    
    
    
[oo,labelsAll]= max(histograms2');
labelsAll = labelsAll>1;
    
    
    
end
    
