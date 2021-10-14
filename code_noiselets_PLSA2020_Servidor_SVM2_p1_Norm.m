%%%%%%%%%%%%
%Using Noiselets to improve Nuclei Segmentation methods
%@version 2.0
%@author M.Sc. Ricardo Moncayo <ramoncayomar@unal.edu.co>
%Ph.D. Eduardo Romero
%%%%%%%%%%%%



addpath('/mnt/md0/ricardo/CodesInternship/GITP/')
addpath('nuclei_seg/staining_normalization/')
addpath('SVM_Noiselets/')
addpath('FNT')
addpath('PLSA_TEM_EM/')
addpath('nuclei_seg')
addpath('nuclei_seg/veta_watershed/')

%%% Dataset Dir %%%
FolderIMG = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Tissue_images/';
FolderXML = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Annotations/';
FolderMasksTest = [FolderIMG,'Masks/'];


FolderTesting = FolderIMG; 
FolderTraining = '/mnt/md0/Histopathology/Datasets/TCIA/'; 
 



%% Folder Results
FolderResults = '//mnt/md0/ricardo/NoiseletProject/Results_SVM_2021/';
mkdir(FolderResults)
FolderImgsMethod = [FolderResults,'ImgsMethodRed2/'];
mkdir(FolderImgsMethod);
FolderImgsOriginal = [FolderResults,'ImgsOriginalRed2/'];
mkdir(FolderImgsOriginal);
FolderImgsSumado  = [FolderResults,'ImgsSumadoRed2/'];
mkdir(FolderImgsSumado);
FolderImgsSinRuido  = [FolderResults,'ImgsSinRuidoRed2/'];
mkdir(FolderImgsSinRuido);
FolderImgsGroundTruth  = [FolderResults,'ImgsGrounTruth2/'];
mkdir(FolderImgsGroundTruth);


%%%%%%%OPTIONS
load([FolderResults,'MatrizExperimentosSVM.mat'])%% archivo que contiene los parametros experimentales
for Expe =1:4%6
% Load parameters for each experiment
Scales = [ExperimentosNoiselets{Expe,2}];
WinPlsa = ExperimentosNoiselets{Expe,3};
ResultsName = ExperimentosNoiselets{Expe,1}{1};
K_clusters =ExperimentosNoiselets{Expe,4};
fprintf('Testing: %s \n',ResultsName)

%%%%%%%Training
ListIMGALL = dir(strcat(FolderTraining,'*_crop.png'));
ListImgsTrain = ListIMGALL;
%%% testing 




%% Feature Extration For dictionary Building
allFeatures = [];
allLabels = [];

VocabCases = randperm(length(ListIMGALL),round(length(ListIMGALL)*0.35));

ListIMGALL = ListIMGALL(VocabCases);

for Lo = 1:length(ListIMGALL)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListIMGALL(Lo).name) 
    ImTest = imread([FolderTraining,ListIMGALL(Lo).name]);
    fprintf('Mejorando Imagen\n') 
 %ToDo: COncatenate img features
    [XFeatures,CoordIMG] = NoiseletsFeatures(ImTest,Scales,WinPlsa);
    allFeatures = [allFeatures;XFeatures];
%    allLabels = [allLabels;imLabels];
    
       
end    



%%% build the dictionary using kmeans

[idx,vocab] = kmeans(allFeatures,K_clusters,'distance','cityblock');
%Build histograms and extract labels
histo_img=[];
labels_img=[];
for Lo = 1:length(ListImgsTrain)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListImgsTrain(Lo).name) 
    fprintf('%d\n',Lo)
    ImTest = imread([FolderTraining,ListImgsTrain(Lo).name]);
    fprintf('Mejorando Imagen\n') 
        %% abrir mascara manual para generar etiquetas
    ImGroundT = imread([FolderTraining,ListImgsTrain(Lo).name(1:end-8),'labeled_mask_corrected.png'])>0;
  %  [binary_maskTest,color_maskTest] =
  %  xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
  %  dataset Monuseg
   [labelsAll,histograms] = NoiseletsPLSAHistogramImg(ImTest,Scales,WinPlsa,ImGroundT,vocab);
       histo_img = [histo_img;histograms];
       labels_img = [labels_img;labelsAll];
       
end



%% train a svm
X=histo_img;
Y=labels_img;

ClassTreeEns = fitensemble(X,Y,'AdaBoostM1',100,'Tree');
ClassModel =ClassTreeEns;
rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%--------- TEST ------- %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Save results


Resultados = cell(length(ListIMGALL),5);
ResultadosSumado =cell(length(ListIMGALL),5);
ResultadosOriginal=cell(length(ListIMGALL),5);

%%
ListImgsTest = dir(strcat(FolderTesting,'*.tif'));

histo_img=[];
labels_imgTest=[];
for Lo = 1:length(ListImgsTest)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListImgsTest(Lo).name) 
    ImTest = imread([FolderTesting,ListImgsTest(Lo).name]);
    fprintf('Mejorando Imagen\n') 
    MaskName = [FolderMasksTest,ListImgsTest(Lo).name];
        %% abrir mascara manual para generar etiquetas
tic
    if ~exist(MaskName)
        fprintf('generating image ')
    [ImGroundT,color_mask]=xlmToMask(ListImgsTest(Lo).name(1:end-4),FolderXML,FolderTesting);
    imwrite(ImGroundT,MaskName)
    
    else
       ImGroundT = imread(MaskName)>0;       
        
    end
    ImManualMask=ImGroundT;
toc
  % [binary_maskTest,color_maskTest] =
  % xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
  % dataset Monuseg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% Removing Noise of the classified histograms%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
%% 
   tic
  try
  [labelsAll,histograms,ImEnhance,ImNormNoiselet] = NoiseletsPLSAHistogramImg(ImTest,Scales,WinPlsa,ImGroundT,vocab,ClassModel);
    catch
        ImEnhance = ImTest;
        fprintf('No se pudo encontrar señal\n')
    end
 
   toc
   %%%%Sacando Mascaras
   fprintf('Sacando Mascara\n') 
   MaskEvaluate = getWatershedMask(ImEnhance);
   MaskOriginal = getWatershedMask(ImTest);
   fprintf('Suma de las mascaras \n')
   sumaBinary = or(MaskOriginal,MaskEvaluate);
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Sin Metodo\n') 
   DiceCoeff = dice(double(MaskEvaluate>0), double(ImManualMask>0)); 
   JAC = jaccard(double(MaskEvaluate>0), double(ImManualMask>0)); 
    Resultados{Lo,1} = ListIMGALL(Lo).name(1:end-4);
    Resultados{Lo,2} = 0;
    Resultados{Lo,3} = DiceCoeff;
    Resultados{Lo,4} = JAC;
    ResultadosOriginal{Lo,5} = 0;
    save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    
   %%%%%%%%%%%%%%%%%%
      fprintf('Calculando Metricas Watershed Original\n') 

   DiceCoeff = dice(double(MaskOriginal>0), double(ImManualMask>0)); 
   JAC = jaccard(double(MaskOriginal>0), double(ImManualMask>0)); 
   ResultadosOriginal{Lo,1} = ListIMGALL(Lo).name(1:end-4);
   ResultadosOriginal{Lo,2} = 0;
   ResultadosOriginal{Lo,3} = DiceCoeff;
   ResultadosOriginal{Lo,4} = JAC;
   ResultadosOriginal{Lo,5} = 0;
   save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    

   
   ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
   imwrite(ImBorde,[FolderImgsOriginal,ListIMGALL(Lo).name(1:end-4),'.png' ]);

   
   
   
   fprintf('Calculando Metricas Suma Metodos\n') 

   
    DiceCoeff = dice(double(sumaBinary>0), double(ImManualMask>0)); 
    JAC = jaccard(double(sumaBinary>0), double(ImManualMask>0)); 
    ResultadosSumado{Lo,1} = ListIMGALL(Lo).name(1:end-4);
    ResultadosSumado{Lo,2} = 0;
    ResultadosSumado{Lo,3} = DiceCoeff;
    ResultadosSumado{Lo,4} = JAC;
    ResultadosSumado{Lo,5} = 0;
    save([FolderResults,'Results_',ResultsName,'Sumados.mat'],'ResultadosSumado')    

    ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));
    imwrite(ImBorde,[FolderImgsSumado,ListIMGALL(Lo).name(1:end-4),'.png' ]);
   
   
   
    %%Imagen borde nucleos detectados en la original
    ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
    imwrite(ImBorde,[FolderImgsOriginal,ListIMGALL(Lo).name(1:end-4),'.png' ]);


    
    %%%%% Imagen sin ruido
  imwrite(ImEnhance,[FolderImgsSinRuido,ListIMGALL(Lo).name(1:end-4),'.png' ]);

    
    
    
  histo_img = [histo_img;histograms];
  labels_imgTest = [labels_imgTest;labelsAll];

       
end


%%%% SACAR RESULTADOS GENERALES DE LA CLASSIFICACION

%% Evaluate SVM
Y_predict=predict(ClassTreeEns,histo_img);

%Results



[c,cm,ind,per] = confusion(labels_imgTest',Y_predict');
accuaracy = sum(diag(cm))/sum(cm(:));
fprintf('accuaracy:%d\n',accuaracy)

save([FolderResults,'Results_',ResultsName,'_SVM.mat'],'accuaracy','cm')    


end



