%%%%%%%%%%%%
%Using Noiselets to improve Nuclei Segementation methods
%@version 1.0
%@author Ricardo Moncayo <ramoncayomar@unal.edu.co>, Christian Arias,
%Eduardo Romero
%%%%%%%%%%%%



nSample = 3;%size of noiselets matrix 2^nSample

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

FolderMASKDenoised = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_mask_Denoised/';
mkdir(FolderMASKDenoised)

FolderMASKGround = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_mask_GroundT/';
mkdir(FolderMASKGround)

FolderMASK = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/Tissue_mask_Baseline/';
mkdir(FolderMASK)

FolderFinalI = '/home/ricardo/Documentos/TESISDOCT/NucleosMTC/final_images/';
mkdir(FolderFinalI)
%%%%%%%Training
ListIMGALL = dir(strcat(FolderIMG,'*.tif'));
Resultados = {};

%%leave one out
for Lo = 1:length(ListIMGALL)
ListIMG=ListIMGALL;
StTest = ListIMG(Lo);
fprintf('%d of %d: %s\n',Lo,length(ListIMGALL),StTest.name);
ListIMG(Lo)=[];

XnucleosAll =[];
XfondoAll = [];
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
        fprintf('Finding Watershed Mask\n');
        M = getWatershedMask( imgNorm );
        imwrite(M, imgPathMASK);
    else
       M = imread(imgPathMASK);
    end
    
    if ~exist(imgPathMASKGround,'file')
        fprinf('computing groun truth mask\n');
    %%groundTruth Mask
       [binary_mask,color_mask,JAI,DiceCoeff] = he_to_binary_mask(ListIMG(n).name(1:end-4),FolderXML,FolderIMG,M);
            imwrite(binary_mask, imgPathMASKGround);
    else
       binary_mask = imread(imgPathMASKGround);
    end
      %ojo cambiar valores tambien en getwateshed
        [~,ImgH,~] = normalizeStaining(imgNorm,220,0.15);
        fprintf('Computing Noiselets features\n')

        [Xnucleos,Xfondo] = FeaturesNoiselets (ImgH,binary_mask,nSample);
        XnucleosAll = [XnucleosAll;Xnucleos];
        XfondoAll = [XfondoAll;Xfondo];
    
    %%%sacar nosiselets agrupar
    %entrenar clasificador
    %sacar imagen sin "ruido"
    %sacar nucleos
    %validar
    
    
    

end
X = real( [XnucleosAll;XfondoAll]);
Y = [ones(length(XnucleosAll),1);zeros(length(XfondoAll),1)];
%%%Train Classifier
fprintf('training clasiffier\n');
SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');
%%% Extract test features
ImgTest = imread(strcat(StTest.folder,'/',StTest.name));
[~,ImgTestH,~] = normalizeStaining(ImgTest,220,0.15);
ImgTestFinal = im2double(rgb2gray(ImgTestH));
fprintf('Noiselets features of test image\n');

[Xtest,imAsoc] = FeaturesNoiseletsTest (ImgTestH,nSample);
%%% evaluate patches
Ypred = predict(SVMModel,real(Xtest));
%%denoise patches
for pre=1:length(Ypred)
   if Ypred(pre)==0
       imAsoc(pre,:)=1;
   end
       
end


Window_size = 2^nSample;
imgY = size(ImgTestFinal,1);
imgX = size(ImgTestFinal,2);
porcionX = imgX/Window_size;
porcionY = imgY/Window_size;

%cuadricular

for Xi =1:porcionX
for Yi=1:porcionY
    ImgTestFinal((Yi-1)*Window_size+1:(Yi)*Window_size,(Xi-1)*Window_size+1:(Xi)*Window_size)=reshape(imAsoc((Xi-1)*floor(porcionX)+Yi,:),Window_size,Window_size);
    fprintf('%d\n',(Xi-1)*floor(porcionX)+Yi)     
end
end
imwrite(ImgTestFinal,strcat(FolderFinalI,StTest.name));


fprintf('COmo vamos hasta aqui')
[binary_maskTest,color_maskTest,JAI,DiceCoeff] = he_to_binary_mask(StTest.name(1:end-4),FolderXML,FolderIMG,ImgTestFinal);

Resultados{Lo,1} = StTest.name;
Resultados{Lo,2} = JAI;
Resultados{Lo,3} = DiceCoeff;

save('secondTryResults.mat','Resultados')    
%% extract nuclei
%% metrics 
end



    %names = ["Nearest Neighbors", "Linear SVM", "RBF SVM","Chi2 SVM","Cosine SVM","GHK SVM" "Gaussian Process",
   %          "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
    %         "Naive Bayes", "QDA","LogisticRegression","GradientBoosting","Voting1vs3","Voting1v2","NewVoting","onevsallVoting"]

         
         %%%%Testing image

%%%metric

%%%mesuares

%%final result
