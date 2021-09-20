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
for Expe = 1:length(ExperimentosNoiselets)
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



%%%% Save results


Resultados = cell(length(ListIMGALL),5);
ResultadosSumado =cell(length(ListIMGALL),5);
ResultadosOriginal=cell(length(ListIMGALL),5);

%% Feature Extration For dictionary Building
allFeatures = [];
allLabels = [];

VocabCases = randperm(length(ListIMGALL),round(length(ListIMGALL)*0.35));

ListIMGALL = ListIMGALL(VocabCases);

for Lo = 1:length(ListIMGALL)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListIMGALL(Lo).name) 
    ImTest = imread([FolderIMG,ListIMGALL(Lo).name]);
    fprintf('Mejorando Imagen\n') 
 %ToDo: COncatenate img features
    [XFeatures,CoordIMG] = NoiseletsFeatures(ImTest,Scales,WinPlsa);
    allFeatures = [allFeatures;XFeatures];
    allLabels = [allLabels;imLabels];
    
       
end    



%%% build the dictionary using kmeans
[idx,vocab] = kmeans(allFeatures,K_clusters);
%Build histograms and extract labels
histo_img=[];
Labels_img=[];
for Lo = 1:length(ListImgsTrain)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListImgsTrain(Lo).name) 
    ImTest = imread([FolderIMG,ListImgsTrain(Lo).name]);
    fprintf('Mejorando Imagen\n') 
        %% abrir mascara manual para generar etiquetas
    ImGroundT = imread([FolderIMG,ListImgsTrain(Lo).name(1:end-8),'labeled_mask_corrected.png'])>0;
  %  [binary_maskTest,color_maskTest] =
  %  xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
  %  dataset Monuseg
   [labelsAll,histograms] = NoiseletsPLSAHistogramImg(ImTest,Scales,WinPlsa,groundT,vocab);
       histo_img = [histo_img;histograms];
       labels_img = [labels_img;labelsAll];
       
end



%% train a svm
histo_img=X;
Y=labels_img;

ClassTreeEns = fitensemble(X,Y,'AdaBoostM1',100,'Tree');

rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%--------- TEST ------- %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ListIMGTest = dir(strcat(FolderTesting,'*tiff'));

histo_img=[];
Labels_img=[];
for Lo = 1:length(ListImgsTrain)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListImgsTrain(Lo).name) 
    ImTest = imread([FolderIMG,ListImgsTrain(Lo).name]);
    fprintf('Mejorando Imagen\n') 
        %% abrir mascara manual para generar etiquetas
    [ImGroundT,color_mask]=xlmToMask(ListIMGTest(Lo).name(1:end-4),folderXML,folderIMG);
  % [binary_maskTest,color_maskTest] =
  % xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
  % dataset Monuseg
 
  [labelsAll,histograms] = NoiseletsPLSAHistogramImg(ImTest,Scales,WinPlsa,groundT,vocab);
  histo_img = [histo_img;histograms];
  labels_img = [labels_img;labelsAll];
 
  
  
 
       
end




%% Evaluate SVM
Y_predict=predict(ClassTreeEns,histo_img);

%Results



[c,cm,ind,per] = confusion(labels_img,Y_predict)

fprintf(c)



end



