function [labelsAll,histograms,im,Inm] = NoiseletsPLSAHistogramImg(im,Scales,WinPlsa,GroundT,vocab,ClassModel,UsedDistance)
%%%%%%%%%%%% Obtain histograms of an image by using and vocabulary and if
%%%%%%%%%%%% classifier is avalaible claffified histograms



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




    win_size=Scales;

    [aa, bb, cc] = size(im);
    andregions = ones(aa,bb);
    [Inorm,H,E] = normalizeStaining(im,im,220,0.06);
    %[H,E] = ColorSepCimalab(im);
    %labim = rgb2lab(H);
    ime = double(H(:,:,1));
    ime=(ime-min(ime(:)))/(max(ime(:))-min(ime(:)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Extract Noiselets from tiles%%%%
    if exist('GroundT')
    [pila_patches_orig,coords,lienzo,pila_patches_gray,labels] = ExtractTilesAndNoiseletsLabels(Scales,ime,aa,bb,GroundT);
    else
        [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(Scales,ime,aa,bb);
    end
        %% cords 
    %%%%%%%%%%%%%%%%%%
    
      %  dataset Monuseg
   % [XFeatures,imLabels] = NoiseletsFeaturesLabels(ime,Scales,WinPlsa,groundT);
   
   abs_pila = real(pila_patches_orig); %Se saca la Magnitud
    angle_pila = imag(pila_patches_orig); %SE SACA EL ANGULO, se le su

    %vector a clusterizacion
    X=zeros(size(pila_patches_orig,1),size(pila_patches_orig,2)*2);
    X(:,1:2:end-1)=abs_pila;
    X(:,2:2:end) = angle_pila; 
    XFeatures = X;

    
        [~,idx] = pdist2(vocab,XFeatures,UsedDistance,'Smallest',1);

  
       
       
       lienzo2 =zeros(size(lienzo));


    
    %% Asigna la ocurrencia al tile
    for L = 1:length(coords)

        aux_point = coords(L,:);
        lienzo(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
            idx(L)*ones(win_size,win_size);
        
        if exist('GroundT')
        lienzo2(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
            labels(L)*ones(win_size,win_size);
        end

    end   
   
    
    
    %%%% PLSA Analysis
    
    Nclusters=size(vocab,1);
    [pila_patches_PLSA,coordsPLSA,lienzoPLSA,histograms] = ExtractTilesOneDPyramidal(WinPlsa,lienzo,aa,bb,Nclusters);
    if exist('GroundT')
        %this can be better to high speed ^<--
    [pila_patches_PLSA2,coordsPLSA2,lienzoPLSA2,histograms2] = ExtractTilesOneDPyramidal(WinPlsa,lienzo2+1,aa,bb,2);
    if any(GroundT(:))
    labelsAll = histograms2(:,2)>round(0.10*WinPlsa^2);
    else
        labelsAll = zeros(size(histograms2(:,1)));
    end
%labelsAll = labelsAll>1;
    
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sacar etiquetas si existen%%% 
    
    



if exist('ClassModel')
           lienzo2 =zeros(size(lienzo));

    LabelsHisto = predict(ClassModel,histograms);
        for L = 1:length(coordsPLSA)
            
        aux_point = coordsPLSA(L,:);

        lienzo2(1+WinPlsa*(aux_point(1)-1):(aux_point(1))*WinPlsa,1+WinPlsa*(aux_point(2)-1):(aux_point(2))*WinPlsa) = ...
            LabelsHisto(L)*ones(WinPlsa,WinPlsa);

        end   
    
    Imbi=lienzo2==0;
    Imbi = bwareaopen(Imbi, 2*WinPlsa*WinPlsa,4);
    idxNoiseRemove = find(Imbi==1);
    imCopy = im;
    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
    imR(idxNoiseRemove) = 255;% double(imR(idxNoiseRemove))/205;
    imG(idxNoiseRemove) = 255;%double(imG(idxNoiseRemove))/114;
    imB(idxNoiseRemove) = 255;%double(imB(idxNoiseRemove))/145;
    im(:,:,1)=imR;
    im(:,:,2)=imG;
    im(:,:,3)=imB;
    [Inm,Hnm,Enm] = normalizeStaining(im,imCopy);
    im=imCopy;
    imR = Inm(:,:,1);
    imG = Inm(:,:,2);
    imB = Inm(:,:,3);
    imR(idxNoiseRemove) = 255;% double(imR(idxNoiseRemove))/205;
    imG(idxNoiseRemove) = 255;%double(imG(idxNoiseRemove))/114;
    imB(idxNoiseRemove) = 255;%double(imB(idxNoiseRemove))/145;
    im(:,:,1)=imR;
    im(:,:,2)=imG;
    im(:,:,3)=imB;
   
    
end
    
    
    
end
    
