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

%%% Dataset Dir %%%
FolderIMG = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Tissue_images/';
FolderXML = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Annotations/';
FolderMasksTest = [FolderIMG,'Masks/'];

load('ListImgsMonuseg.mat')

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%OPTIONS
      load([FolderResults,'MatrizExperimentosSVM.mat'])%% archivo que contiene los parametros experimentales
      
      
for ClassExp = 1:length(clasesMonu) %realizar experimento por clase
SelClass = clasesMonu(ClassExp)    ;
fprintf('Evaluating for %s tissue\n',SelClass)
count=0;
    for rn = 1:length(ListOfImgCases)%take out list per tissues
       ItClass = ListOfImgCases{rn,2}==SelClass; % logical existence of a determine class
       
       if ItClass
          count=count+1; 
          imageForExperiment(count) = ListOfImgCases{rn,2} ;
       end
      
    end

           
    %%%% Save Results
    Resultados_acc_SVM = cell(length(imageForExperiment),length(ExperimentosNoiselets));
    Resultados_SumadoImg = cell(length(imageForExperiment),length(ExperimentosNoiselets));
    Resultados_OriginalImg = cell(length(imageForExperiment),length(ExperimentosNoiselets));
    Resultados_Only_method = cell(length(imageForExperiment),length(ExperimentosNoiselets));

    %%%%%%%%%%%%%
    
    
    %%%%% LEAVE ONE OUT VALIDATION SCHEME %%%%
    
    for OneOut = 1:length(imageForExperiment)
       TestImg =  imageForExperiment(OneOut);
       idxTraining = not(imageForExperiment==TestImg);
       ImgsTraining = imageForExperiment(idxTraining);
       
       %%%%%%% Training %%%%%%%%%
       ListIMGALL = ImgsTraining; %just a change of variable
       ListImgsTrain = ListIMGALL;



    for Expe = 13:16%6
        % Load parameters for each experiment
        Scales = [ExperimentosNoiselets{Expe,2}];
        WinPlsa = ExperimentosNoiselets{Expe,3};
        ResultsName = ExperimentosNoiselets{Expe,1}{1};
        K_clusters =ExperimentosNoiselets{Expe,4};
        fprintf('Testing: %s \n',ResultsName)


    %%% testing 

    %% Feature Extration For dictionary Building
    allFeatures = [];
    allLabels = [];

    %VocabCases = randperm(length(ListIMGALL),round(length(ListIMGALL)*0.35));
    VocabCases = randperm(length(ListIMGALL),length(ListIMGALL));

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
  %  ImGroundT = imread([FolderTraining,ListImgsTrain(Lo).name(1:end-8),'labeled_mask_corrected.png'])>0;
    [binary_maskTest,color_maskTest] = xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
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


histo_img=[];
labels_imgTest=[];
for Lo = 1:length(TestImg) %%% this for is always 1, a cicle for her would be used to K-fold validation
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',TestImg(Lo).name) 
    ImTest = imread([FolderTesting,TestImg(Lo).name]);
    fprintf('Mejorando Imagen\n') 
    MaskName = [FolderMasksTest,TestImg(Lo).name];
        %% abrir mascara manual para generar etiquetas
tic
    if ~exist(MaskName)
        fprintf('generating image ')
    [ImGroundT,color_mask]=xlmToMask(TestImg(Lo).name(1:end-4),FolderXML,FolderTesting);
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
   fprintf('Calculando Metricas Sin Sumar\n') 
   DiceCoeff = dice(double(MaskEvaluate>0), double(ImManualMask>0)); 
   Resultados_Only_method{Lo,1} = ListIMGALL(Lo).name(1:end-4);    
   Resultados_Only_method{Lo,2} = DiceCoeff;


   %%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Watershed Original\n') 
   DiceCoeff = dice(double(MaskOriginal>0), double(ImManualMask>0)); 
   Resultados_OriginalImg{Lo,1} = ListIMGALL(Lo).name(1:end-4);    
   Resultados_OriginalImg{Lo,2} = DiceCoeff;


   %save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    

   
 %  ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
  % imwrite(ImBorde,[FolderImgsOriginal,ListIMGALL(Lo).name(1:end-4),'.png' ]);

   
   
   
   fprintf('Calculando Metricas Suma Metodos\n') 

   
    DiceCoeff = dice(double(sumaBinary>0), double(ImManualMask>0)); 
    Resultados_SumadoImg{Lo,1} = ListIMGALL(Lo).name(1:end-4);
    Resultados_SumadoImg{Lo,2} = DiceCoeff;
    
    save([FolderResults,'Results_',SelClass,'_Space.mat'],...,
        'Resultados_SumadoImg','Resultados_acc_SVM','Resultados_Only_method','Resultados_OriginalImg')

    %ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));    
    %imwrite(ImBorde,[FolderImgsSumado,ListIMGALL(Lo).name(1:end-4),'.png' ]);
   
   
   
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
    end
end






