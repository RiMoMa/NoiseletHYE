function [im,im2] = NoiseletsPLSAHistogramImg(im,Scales,WinPlsa,K_clusters,groundT,vocab)
%%% Sacar Mascaras
addpath('PLSA_TEM_EM')
addpath('FNT')
%%% OPTIONS %%%
cluster='kmeans'; % 'plsa'
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





[aa, bb, cc] = size(im);
andregions = ones(aa,bb);
t = 0;
for win_size = Scales
    t = 2; % t=t+1; this is for a normalization without noise and then a noise remotion
    [Inorm,H,E] = normalizeStaining(im,im,220,0.06);
    %[H,E] = ColorSepCimalab(im);
    %labim = rgb2lab(H);
    ime = H(:,:,1);

    %%% Extract Noiselets from tiles
    [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(win_size,ime,aa,bb);
    
    
    lienzo2 = lienzo;

    abs_pila = real(pila_patches_orig); %Se saca la Magnitud
    angle_pila = imag(pila_patches_orig); %SE SACA EL ANGULO, se le su

    %vector a clusterizacion
    X=zeros(size(pila_patches_orig,1),size(pila_patches_orig,2)*2);
    X(:,1:2:end-1)=abs_pila;
    X(:,2:2:end) = angle_pila;
    
    idx=vocab;%TOdo PDIST
    
    %% Asigna la ocurrencia al tile
    for L = 1:length(coords)

        aux_point = coords(L,:);
        lienzo(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
            idx(L)*ones(win_size,win_size);

    end   
   
    
    
    %%%% PLSA Analysis
    
    
    [pila_patches_PLSA,coordsPLSA,lienzoPLSA,histograms] = ExtractTilesOneDPyramidal(WinPlsa,lienzo,aa,bb);
    
    [idx, prob_topic_doc,lls] = plsa(normc(histograms)', 2);
    %[prob_term_topic, prob_topic_doc,lls]
    [s,id2]=max(prob_topic_doc);
    
     for L = 1:length(coordsPLSA)
         aux_point = coordsPLSA(L,:);
        lienzoPLSA(1+WinPlsa*(aux_point(1)-1):(aux_point(1))*WinPlsa,1+WinPlsa*(aux_point(2)-1):(aux_point(2))*WinPlsa) = ...
            id2(L)*ones(WinPlsa,WinPlsa);

    end
    
    

%%%%%%%Encuentra el cluster de ruido medio y señal%%%%%%%%%%%%

    IntensidadTotal = zeros(1,K_clusters);
    % buscar las intensidades de los clusters
    EnH = normc(histograms);
    for bigL =1:2
     IntensidadTotal(1,bigL)  = entropy(EnH(id2==bigL,:));
    end
    [ value,posRuido] = sort(IntensidadTotal,'descend');
    if t==1
    Imbi=lienzoPLSA==posRuido(2);
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
    im = Inm;
elseif t==2
%     BIsignal = lienzoPLSA==posRuido(1);

%     imR = im(:,:,1);
%     imG = im(:,:,2);
%     imB = im(:,:,3);
%     imR(idxSignalEnhance) = imR(idxSignalEnhance)-30;
%     imG(idxSignalEnhance) = imG(idxSignalEnhance)-30 ;
%     imB(idxSignalEnhance) = imB(idxSignalEnhance)-30;
%     im(:,:,1)=imR;
%     im(:,:,2)=imG;
%     im(:,:,3)=imB;
    
    %% Ahora con las Noiselets
    IntensidadTotal2 = zeros(1,K_clusters);
    % buscar las intensidades de los clusters
    for bigL =1:2
     %idx_C = find(lienzo==bigL);%
     W  = imhist((ime.*(uint8(lienzo==bigL))+255*uint8(~(lienzo==bigL))).*uint8(lienzoPLSA==posRuido(1)));
     W = W(10:150);
     IntensidadTotal2(1,bigL) = entropy(W/max(W(:)));
    end
    [ value,posSignalNoiselet] = sort(IntensidadTotal2,'descend');
    BIsignal = lienzo==posSignalNoiselet(1);
    BIsignal = bwareaopen(BIsignal, 2*WinPlsa*4*2,4);
   % idxSignalEnhance = find(BIsignal==1);
    %elEs = strel('disk',1,6);
    %BIsignal = imdilate(BIsignal,elEs);
    idxSignalEnhance = find(BIsignal==1);
    im2  = im;

    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
    imR(idxSignalEnhance) = double(imR(idxSignalEnhance))/255*145;
    imG(idxSignalEnhance) = double(imG(idxSignalEnhance))/255*124 ;
    imB(idxSignalEnhance) = double(imB(idxSignalEnhance))/255*169;
    im(:,:,1)=imR;
    im(:,:,2)=imG;
    im(:,:,3)=imB;
    elEs = strel('disk',4,6);
    BIsignal = imdilate(BIsignal,elEs);
    
    BIsignal2 = ~BIsignal;
    BIsignal2 =BIsignal2- bwareaopen(BIsignal2, 2*WinPlsa*WinPlsa,4);
    BIsignal = or(BIsignal2,BIsignal);
    

    idxSignalEnhance = find(BIsignal==0);
    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
    imR(idxSignalEnhance) = 255;
    imG(idxSignalEnhance) = 255 ;
    imB(idxSignalEnhance) = 255;
    im(:,:,1)=imR;
    im(:,:,2)=imG;
    im(:,:,3)=imB;
    
    imR = im2(:,:,1);
    imG = im2(:,:,2);
    imB = im2(:,:,3);
    imR(idxSignalEnhance) = 255;
    imG(idxSignalEnhance) = 255 ;
    imB(idxSignalEnhance) = 255;
    im2(:,:,1)=imR;
    im2(:,:,2)=imG;
    im2(:,:,3)=imB;
    
    
    
end
    
end


end
    % ordenar las intensidades
%     [IntensidadTotal,sortKdetermined] = sort(IntensidadTotal,'descend');
%     
% %%%%%%%%%%%%%%%    %%%%%%%%%%%%%%
%    
%     imeWithoutNoise = ime;%copia de la imagen analizada
%     %%%%%%%%%%%% Recorre Parches Para eliminar ruido
%     for Nr = 1:size(IntensidadTotal,2)
%      % eliminar ruido grande del cluster de ruido
%         lienzoN = lienzo==sortKdetermined(Nr); %encuentra los parches del cluster Nr
%         regions = bwlabel(lienzoN); % genera una imagen regiones
%         areas = regionprops(regions,'area','EquivDiameter'); % extrae el area encontrada de cada region
%         idxAreas = [1:length(areas)]; % Cantidad de areas
%         accepted = ruleBasedClassification(areas,{'Area', [win_size^2-1 60*60],'EquivDiameter',[1, 35]} );% eliminar areas grandes y pequeñas
%         idx_rejected = idxAreas(~accepted); % cuales son las areas rechazadas
%         idx_accepted = idxAreas(accepted);
%         if Nr==2
%             idx_accepted = idxAreas;
%             idx_rejected = [];
%         end
%         % Aqui se relaciona las coords con el numero de region%%
%         regionsOfCords = zeros(1,size(coords,1));
%         for L = 1:length(coords)
%             aux_point = coords(L,:);
%             regionsOfCords(L) = unique(regions(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) );
%         end    
%     %%%%%%
%     
%     %%% Aca se va a eliminar las regiones con areas grandes de la imagen
%     %%% binaria o de regiones
%     
%     lienzoNCopy = lienzoN;
%     for re = 1:length(idx_rejected)
%         ToRemove = regions==idx_rejected(re); %encuentra la ROI a eliminar
%         lienzoNCopy = lienzoNCopy-ToRemove; % la elimina de la imagen de regiones del cluster Nr
%     end
%     CoordsAcepted =[];
%     for re = 1:length(idx_accepted)
%         CorRegion = find(regionsOfCords==idx_accepted(re));%encuentra que Coords fue eliminado
%         CoordsAcepted =[CoordsAcepted,CorRegion]; %concatena todos los coords
%     end
%     %%% Remueve el ruido o el area grande del cluster 1 y 2
%     NoiseImage = lienzoN - lienzoNCopy;
%     RemainImage = lienzoNCopy;
%     idx_toremove = find(NoiseImage==1);
%     imeWithoutNoise(idx_toremove)=255;
%     
%     %%% Busco las noiselets que aun quedan para volver a clusterizar
%     BinaryCords = logical(ones(1,size(coords,1)));
%     
%         BinaryCords(CoordsAcepted)=0;
%         BinaryCords = ~BinaryCords;
%     
%     Xremain = X(BinaryCords,:);
%     idxRemain = kmeans(Xremain,2,'distance','cityblock');
%     CoordsRemain = coords(BinaryCords,:);
%     lienzo2 = zeros(size(lienzo));
%    for L = 1:length(CoordsRemain)
%        aux_point = CoordsRemain(L,:);
%        lienzo2(1+win_size*(aux_point(1)-1):(aux_point(1))*win_size,1+win_size*(aux_point(2)-1):(aux_point(2))*win_size) = ...
%             idxRemain(L)*ones(win_size,win_size);
%    end
%     
%     if Nr ==1
%         RemainImageC1 = RemainImage;
%         
%         % en solo ruido conservo el cluster que tenga mas areas grandes (info de nucleo)
%         perimeter1 = median(struct2array(regionprops(lienzo2==1,'perimeter')));
%         perimeter2 = median(struct2array(regionprops(lienzo2==2,'perimeter')));
%       
%     if perimeter1>perimeter2
%            lienzoNCopy = lienzo2==1;
%     else
%         lienzoNCopy = lienzo2==2;
%     
%     end
%     
%         %%% Remueve el ruido del cluster determinado y genera una imagen
%         %%% total de ruido
%         %NoiseImage = lienzoNCopy;
%         
%         NoiseImage2 = NoiseImage+(RemainImage - lienzoNCopy);
%     
%     elseif Nr==2
%         RemainImageC2 = RemainImage;
% %%%%%%%Encuentra el cluster de ruido medio y señal%%%%%%%%%%%%
%    IntensidadTotalCluster2 = zeros(1,2);
%     % buscar las intensidades de los clusters
%     for bigL =1:2
%       idx_C = find(lienzo2==bigL);%
%       IntensidadTotalCluster2(1,bigL)  = entropy(ime(idx_C));
% 
%     end
%     
%     [ value,posNr2] = max(IntensidadTotalCluster2 );
%     RemainImage2 = lienzo2 == posNr2;
%     elseif Nr==3
%         RemainImageC3 = RemainImage;
%         RemainImage3 = lienzoN;
%     end
%     
%    % region2 = regions==idx_area(1);
%    % andregions = and(andregions,lienzoN);
% end
% 
% %idx_toremove = find(NoiseImage2==1);
% %RemainTotal = RemainImage2+RemainImage3;
% %figure;imshow(imoverlay(im,RemainTotal,'cyan'),'InitialMagnification',67)
% idxSignalEnhance = find(RemainImage3==1);
% idxSignalEnhance2 = find(RemainImage2==1);
% 
% idxNoiseRemove = find(NoiseImage2==1);
% imCopy = im;
% imR = im(:,:,1);
% imG = im(:,:,2);
% imB = im(:,:,3);
% imR(idxNoiseRemove) = 255;
% imG(idxNoiseRemove) = 255;
% imB(idxNoiseRemove) = 255;
% im(:,:,1)=imR;
% im(:,:,2)=imG;
% im(:,:,3)=imB;
% imR(idxSignalEnhance) = imR(idxSignalEnhance)-30;
% imG(idxSignalEnhance) = imG(idxSignalEnhance)-30 ;
% imB(idxSignalEnhance) = imB(idxSignalEnhance)-30;
% im(:,:,1)=imR;
% im(:,:,2)=imG;
% im(:,:,3)=imB;
% imR(idxSignalEnhance2) = imR(idxSignalEnhance2)-60;
% imG(idxSignalEnhance2) = imG(idxSignalEnhance2)-60 ;
% imB(idxSignalEnhance2) = imB(idxSignalEnhance2)-60;
% im(:,:,1)=imR;
% im(:,:,2)=imG;
% im(:,:,3)=imB;
% 
% %imeWithoutNoise = im(:,:,1);
% %imeWithoutNoise(idx_toremove)=255;
% %imeWithoutNoiseBlurred = impyramid(impyramid(imeWithoutNoise,'reduce'),'expand');
% %MwithoutNoiseBlurred = getWatershedMask2Noiselet(im,imeWithoutNoiseBlurred);
% %Morigina = getWatershedMask(im);
% 
% %BW = boundarymask(MwithoutNoiseBlurred);
% %figure;imshow(imoverlay(im,BW,'cyan'),'InitialMagnification',67)
% %BW = boundarymask(Morigina);
% %figure;imshow(imoverlay(im,BW,'cyan'),'InitialMagnification',67)
% %Mout = or(MwithoutNoiseBlurred,Morigina);
% %BW = boundarymask(Mout);
% %figure;imshow(imoverlay(im,BW,'cyan'),'InitialMagnification',67)
% 
% 
% end