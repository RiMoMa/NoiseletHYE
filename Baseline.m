%%% CHallengue to England
%Baseline

cd /home/ricardo/Documentos/TESISDOCT/NucleosMTC
addpath(genpath('/home/ricardo/Documentos/TESISDOCT/NucleosMTC/codes/regist_dist_color'))
addpath(genpath('/home/ricardo/Documentos/TESISDOCT/NucleosMTC/codes/nuclei_seg'))

FolderIMG = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_images/';
FolderXML = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Annotations/';


%%%OPEN REFERENCE IMAGE For Normalizations
imgtoNorm = imread('/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_images/TCGA-AR-A1AS-01Z-00-DX1.tif');
%imgtoNorm = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_images/TCGA-HE-7129-01Z-00-DX1.tif' ;
FolderIMGNORM = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_images_NORM/';
mkdir(FolderIMGNORM)

FolderMASK = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_mask_Baseline/';
mkdir(FolderMASK)

FolderMASKGround = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_mask_GroundT/';
mkdir(FolderMASKGround)


ListIMG = dir(strcat(FolderIMG,'*.tif'));
Resultados = {};

for n = 1:length(ListIMG)
    imgPath = strcat(ListIMG(n).folder,'/',ListIMG(n).name);
    imgPathNorm = strcat(FolderIMGNORM,ListIMG(n).name);  
    imgPathMASK = strcat(FolderMASK,ListIMG(n).name);
    imgPathMASKGround = strcat(FolderMASKGround,ListIMG(n).name);
 
    
    %Normalization
    
    if ~exist(imgPathNorm,'file')
        img = imread(imgPath);
        imgNorm = normalizeCPRstain(img, imgtoNorm);
        imwrite(imgNorm, imgPathNorm);
    else
        imgNorm = imread(imgPathNorm);
        
    end    
    
    %Color Deconvolution%nuclei Detection
        
    if ~exist(imgPathMASK,'file')
        M = getWatershedMask( imgNorm );
        imwrite(M, imgPathMASK);
    else
       M = imread(imgPathMASK);
    end
        %%groundTruth Mask
       [binary_mask,color_mask,JAI,DiceCoeff] = he_to_binary_mask(ListIMG(n).name(1:end-4),FolderXML,FolderIMG,M);
    
Resultados{n,1} = ListIMG(n).name;
Resultados{n,2} = JAI;
Resultados{n,3} = DiceCoeff;

save('BaselineResults.mat','Resultados')    

end