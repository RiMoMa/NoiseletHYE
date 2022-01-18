function [labelsAll,histograms,im,Inm] = NoiseletsPLSAHistogramImg_overlap(im,Scales,WinPlsa,GroundT,vocab,ClassModel)
%%%%%%%%%%%% Obtain histograms of an image by using and vocabulary and if
%%%%%%%%%%%% classifier is avalaible claffified histograms

%%% Outputs %%%
Join_labesAll = cell(4,1);
Join_histograms = cell(4,1);
Join_lienzos = cell(4,1);
Join_im = cell(4,1);
Join_Inm = cell(4,1);
%%%%


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
    Control_im = im;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Extract Noiselets from tiles%%%%
    
    for shift = 1:4
        if shift == 1
            % original
            ime_c = ime;
            
            
        elseif shift == 2
            %down shift
            ime_c = ime(win_size/2+1:end,:);
            im=cat(3,im(win_size/2+1:end,:,1),im(win_size/2+1:end,:,2),im(win_size/2+1:end,:,3));
        elseif shift == 3
            %right shift
            ime_c = ime(:,win_size/2+1:end);
            im = cat(3,im(:,win_size/2+1:end,1),im(:,win_size/2+1:end,2),im(:,win_size/2+1:end,3));
            
        elseif shift == 4
            ime_c= ime(win_size/2+1:end,win_size/2+1:end);         
            im = cat(3,im(win_size/2+1:end,win_size/2+1:end,1),im(win_size/2+1:end,win_size/2+1:end,2),im(win_size/2+1:end,win_size/2+1:end,3));
            
        else
            fprintf ('something goes too bad\n')
        end        

    if exist('GroundT') %% if is need to obtain labels for train
        if shift == 1
            % original
            GroundT_c = GroundT;
            
        elseif shift == 2
            %down shift
            GroundT_c = GroundT(win_size/2+1:end,:);
        elseif shift == 3
            %right shift
            GroundT_c = GroundT(:,win_size/2+1:end);
            
        elseif shift == 4
            GroundT_c= GroundT(win_size/2+1:end,win_size/2+1:end);          
            
        else
            fprintf ('something goes too bad\n')
        end
               
        
        [aa, bb, cc] = size(ime_c);
        [pila_patches_orig,coords,lienzo,pila_patches_gray,labels] = ExtractTilesAndNoiseletsLabels(Scales,ime_c,aa,bb,GroundT_c);
       
    else % or only features for test
        [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(Scales,ime_c,aa,bb);
       
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
    [~,idx] = pdist2(vocab,XFeatures,'cityblock','Smallest',1);
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
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sacar etiquetas si existen%%% 
    
    
    %%%% put tiles
    
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
    



    if exist('ClassModel')
        lienzo2 =zeros(size(lienzo));
        LabelsHisto = predict(ClassModel,histograms); % TODO: Probar con probabilidades
        for L = 1:length(coordsPLSA)
            aux_point = coordsPLSA(L,:);
            lienzo2(1+WinPlsa*(aux_point(1)-1):(aux_point(1))*WinPlsa,1+WinPlsa*(aux_point(2)-1):(aux_point(2))*WinPlsa) = ...
                LabelsHisto(L)*ones(WinPlsa,WinPlsa);
        end
        if shift ==1
          %  Imbi=lienzo2==0;            
            Imbi=lienzo2;
            
        elseif shift == 2            
            %Imbi2=lienzo2==0;
            Imbi2=lienzo2;
            lienzo_con(win_size/2+1:end,:) = Imbi2;
            Imbi = or(Imbi,lienzo_con);            
            
        elseif shift == 3                        
            Imbi2=lienzo2;
            lienzo_con(:,win_size/2+1:end) = Imbi2;
            Imbi = or(Imbi,lienzo_con);
        elseif shift ==4
            Imbi2=lienzo2;
            lienzo_con(win_size/2+1:end,win_size/2+1:end) = Imbi2;
            Imbi = or(Imbi,lienzo_con);
            Imbi=Imbi==0;
            %  Imbi = bwareaopen(Imbi, 2*WinPlsa*WinPlsa,4);
        idxNoiseRemove = find(Imbi==1);
        imCopy = im;
        im = Control_im;
        imR = im(:,:,1);
        imG = im(:,:,2);
        imB = im(:,:,3);
        imR(idxNoiseRemove) = 255;% double(imR(idxNoiseRemove))/205;
        imG(idxNoiseRemove) = 255;%double(imG(idxNoiseRemove))/114;
        imB(idxNoiseRemove) = 255;%double(imB(idxNoiseRemove))/145;
        im(:,:,1)=imR;
        im(:,:,2)=imG;
        im(:,:,3)=imB;
        [Inm,Hnm,Enm] = normalizeStaining(im,Control_im);
        im=Control_im;
        imR = Inm(:,:,1);
        imG = Inm(:,:,2);
        imB = Inm(:,:,3);
        imR(idxNoiseRemove) = 255;% double(imR(idxNoiseRemove))/205;
        imG(idxNoiseRemove) = 255;%double(imG(idxNoiseRemove))/114;
        imB(idxNoiseRemove) = 255;%double(imB(idxNoiseRemove))/145;
        im(:,:,1)=imR;
        im(:,:,2)=imG;
        im(:,:,3)=imB;
        Join_im{shift} = im;
        Join_Inm{shift} = Inm;
        
        end
        lienzo_con = zeros(size(Imbi));
      
        
            
    end
    Join_labesAll{shift} = labelsAll;
    Join_histograms{shift} = histograms;

    im = Control_im;
    
    end
    
    labelsAll = cell2mat(Join_labesAll);
    histograms = cell2mat(Join_histograms);
    im = Join_im{shift};
    Inm = Join_Inm{shift};
    
end
    
