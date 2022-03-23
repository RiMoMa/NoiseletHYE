%%%%%%%%%%%%
%Using Noiselets to improve Nuclei Segmentation methods
%@version 2.0
%@author M.Sc. Ricardo Moncayo <ramoncayomar@unal.edu.co>
%Ph.D. Eduardo Romero
%%%%%%%%%%%%


%%%% Libraries %%%%%%%%%%%5
addpath('/mnt/storage/ramoncayomar/Codes/CodesInternship/GITP/')
%addpath('/home/ricardo/Documents/Doctorado/SYNC/all_codes/CodesInternship/GITP/')

addpath('nuclei_seg/staining_normalization/')
addpath('SVM_Noiselets/')
addpath('FNT')
addpath('PLSA_TEM_EM/')
addpath('nuclei_seg')
addpath('nuclei_seg/veta_watershed/')
%%%%%%%%%%%%%%%%%5

%%%%%%%%% Folder Results %%%%%%%%
%FolderResults =     '/home/ricardo/Documents/Doctorado/NucleosMTC/ResultsTCIA_Norm_2022/'; %folder Out
FolderResults =     '/mnt/storage/ramoncayomar/ResultsNoiselets_Original_Complete_noOverlap_2602/'; %folder Out
mkdir(FolderResults)
copyfile('MatrizExperimentosSVM_2022.mat',FolderResults)
FolderImgsMethod = [FolderResults,'ImgsMethodRed2/']; % imgs H&E without noise
mkdir(FolderImgsMethod);
FolderImgsOriginal = [FolderResults,'ImgsOriginalRed2/']; %% make a copy of the data
mkdir(FolderImgsOriginal);
FolderImgsSumado  = [FolderResults,'ImgsSumadoRed2/']; %% borderline of the detected nuclei when mask are added
mkdir(FolderImgsSumado);
FolderImgsSinRuido  = [FolderResults,'ImgsSinRuidoRed2/'];%%% Imgs without Noise masks
mkdir(FolderImgsSinRuido);
FolderImgsGroundTruth  = [FolderResults,'ImgsGrounTruth2/']; %% borderline groundtruth
mkdir(FolderImgsGroundTruth);
%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Dataset Dir %%%
%FolderIMG = '/mnt/md0/Histopathology/Datasets/TCIA/';
FolderIMG = '/mnt/storage/ramoncayomar/Datasets/TCIA/';   
%FolderIMG = '/home/ricardo/Documents/Doctorado/DatasetTCIA_NORM/';
%FolderXML =
%'/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Annotations/'; %this is for monuseg dataset
FolderMasksTest = [FolderIMG,'Masks/'];
FolderTraining = FolderIMG;
FolderTesting = FolderIMG; %change for future datasets
load('TCIADetail.mat')
ListOfImgCases = DetailData;

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%OPTIONS
load([FolderResults,'MatrizExperimentosSVM_2022.mat'])%% archivo que contiene los parametros experimentales
      
      
for ClassExp = 1%:length(clasesMonu) %realizar experimento por clase
    imageForExperiment=[];

    SelClass = clasesMonu(ClassExp) ;
    fprintf('Evaluating for %s tissue\n',SelClass)
    count=0; 
    
    for rn = 1:length(ListOfImgCases)%take out list per tissue selection
       ItClass = ListOfImgCases{rn,2}==SelClass; % logical existence of a determine class
       
       if ItClass
          count=count+1; 
          imageForExperiment(count) = ListOfImgCases{rn,1} ;
       end
    end

           
    %%%% Save Results
    Resultados_acc_SVM = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_SumadoImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_OriginalImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_Only_method = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    
    %%%% PRECISION DETECTION PER CASE
    Resultados_PrecisionD_SumadoImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_PrecisionD_OriginalImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_PrecisionD_Only_method = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    %%%% RECALL DETECTION PER CASE    
    Resultados_RecallD_SumadoImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_RecallD_OriginalImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_RecallD_Only_method = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    %%%% FSCORE DETECTION PER CASE    
    Resultados_FscoreD_SumadoImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_FscoreD_OriginalImg = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    Resultados_FscoreD_Only_method = cell(length(ExperimentosNoiselets),length(imageForExperiment)*2);
    
    
    %%%%%%%%%%%%% K-Fold Partition
    %%% 1. Randomize Images
    %%% 2. Determine number of train and test Images
    KFold = 10;
    cv = cvpartition(imageForExperiment,'KFold',KFold);
    %cv = cvpartition(imageForExperiment,'Holdout',0.3)
    %%% 3. Make data partion
    
    
    %%%%% K FOLD VALIDATION SCHEME %%%%
    
    AccumTest=0; %This variable is employed to count evaluated images and then create the matrix of results dependend of the kfold and the number of test images   
    
    for OneOut = 1:KFold
        
       TestImg = imageForExperiment(cv.test(OneOut));
       
       idxTraining = cv.training(OneOut);
       ImgsTraining = imageForExperiment(idxTraining);
       
       %%%%%%% Training %%%%%%%%%
       ListIMGALLT = ImgsTraining; %just a change of variable
       ListImgsTrain = ListIMGALLT;



    for Expe = 1:length(ExperimentosNoiselets)%16%6
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

    VocabCases = randperm(length(ListIMGALLT),round(length(ListIMGALLT)*0.4));
    %VocabCases = randperm(length(ListIMGALL),length(ListIMGALL));

    ListIMGALL = ListIMGALLT(VocabCases);
%%% VOCABULARY FEATURE EXTRACTION
for Lo = 1:length(ListIMGALL)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',num2str(ListIMGALL(Lo))) 
    ImTest = imread([FolderTraining,num2str(ListIMGALL(Lo)),'_crop.png']);
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
    ImgTrainName = num2str(ListImgsTrain(Lo));
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ImgTrainName) 
    
    ImTrain = imread([FolderTraining,ImgTrainName,'_crop.png']);
    fprintf('Mejorando Imagen\n') 
  %% abrir mascara manual para generar etiquetas
   ImGroundT = imread([FolderTraining,ImgTrainName,'_labeled_mask_corrected.png'])>0;
%    [binary_maskTest,color_maskTest] = xlmToMask(ImgTrainName,FolderXML,FolderIMG);%% Para
  %  dataset Monuseg
    % ImGroundT = binary_maskTest>0;
      [labelsAll,histograms] = NoiseletsPLSAHistogramImg(ImTrain,Scales,WinPlsa,ImGroundT,vocab);
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

fprintf('Test Img\n')
histo_img=[];
labels_imgTest=[];
for Lo = 1:length(TestImg) %%% this for is always 1, a cicle for her would be used to K-fold validation
    AccumTest = AccumTest+1;
    fprintf('AbriendoImagen\n') 
    ImgTestS = num2str(TestImg(Lo));
    fprintf('%s\n',ImgTestS) 
    ImTest = imread([FolderTesting,ImgTestS,'_crop.png']);
    fprintf('Mejorando Imagen\n') 
    MaskName = [FolderMasksTest,ImgTestS,'_labeled_mask_corrected.png'];
        %% abrir mascara manual para generar etiquetas
tic
    if ~exist(MaskName)
        fprintf('generating image \n')
    ImGroundTE = imread([FolderTraining,ImgTestS,'_labeled_mask_corrected.png']);
    ImGroundT = ImGroundTE >0; 
        %[ImGroundT,color_mask]=xlmToMask(ImgTestS,FolderXML,FolderTesting);
    
  %  imwrite(ImGroundT,MaskName)
    
    else
       ImGroundTE = imread(MaskName);       
       ImGroundT = ImGroundTE >0; 
    end
    ImManualMask=ImGroundT;
toc
  % [binary_maskTest,color_maskTest] =
  % xlmToMask(ImgTestS,FolderXML,FolderIMG);%% Para
  % dataset Monuseg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% Removing Noise of the classified histograms%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
%% 
   tic
  try
  [labelsAll,histograms,ImEnhance,ImNormNoiselet] = NoiseletsPLSAHistogramImg(ImTest,Scales,WinPlsa,ImGroundT,vocab,ClassModel);
         %%%% SACAR RESULTADOS GENERALES DE LA CLASSIFICACION

        %% Evaluate SVM
        Y_predict=predict(ClassTreeEns,histograms);

        %Results



        [c,cm,ind,per] = confusion(labelsAll',Y_predict');
        accuaracy = sum(diag(cm))/sum(cm(:));
        fprintf('accuaracy:%d\n',accuaracy)

    
    catch
        ImEnhance = ImTest;
        accuaracy = 0; % this value was not computed
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
%%%%%%%%%%%%%%%%%%%%%%%%% METRICS %%%%%%%%%%%
   fprintf('Calculando Metricas Sin Sumar\n') 
   DiceCoeff = dice(double(MaskEvaluate>0), double(ImManualMask>0)); 
   %% DETECTION METRIC
   Resultados_Only_method{Expe,2*AccumTest-1} = ImgTestS;    
   Resultados_Only_method{Expe,AccumTest*2} = DiceCoeff;
  

   %%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Watershed Original\n') 
   DiceCoeff = dice(double(MaskOriginal>0), double(ImManualMask>0)); 
   Resultados_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;    
   Resultados_OriginalImg{Expe,AccumTest*2} = DiceCoeff;

   %save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    
   
 %  ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
  % imwrite(ImBorde,[FolderImgsOriginal,ImgTestS,'.png' ]);
  
   
   fprintf('Calculando Metricas Suma Metodos\n') 
  
    DiceCoeff = dice(double(sumaBinary>0), double(ImManualMask>0)); 
    Resultados_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_SumadoImg{Expe,AccumTest*2} = DiceCoeff;
 
    
  %  save([FolderResults,'Results_',ResultsName,'_SVM.mat'],'accuaracy','cm')    


    Resultados_acc_SVM{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_acc_SVM{Expe,AccumTest*2} = accuaracy;
    
    
    %%% Metrics for nuclei Detection
 
[Precision_S,Recall_S,Fscore_S] = Nuclei_centroidsBased_metric(ImGroundTE, sumaBinary);
[Precision_O,Recall_O,Fscore_O] = Nuclei_centroidsBased_metric(ImGroundTE, MaskOriginal);
[Precision_M,Recall_M,Fscore_M] = Nuclei_centroidsBased_metric(ImGroundTE, MaskEvaluate);

       %%%% PRECISION DETECTION PER CASE
    Resultados_PrecisionD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_SumadoImg{Expe,AccumTest*2} = Precision_S(3);
    
    Resultados_PrecisionD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_OriginalImg{Expe,AccumTest*2} = Precision_O(3);
    
    Resultados_PrecisionD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_Only_method{Expe,AccumTest*2} = Precision_M(3);
        
    %%%% RECALL DETECTION PER CASE    
    Resultados_RecallD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_SumadoImg{Expe,AccumTest*2} = Recall_S(3);
    
    Resultados_RecallD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_OriginalImg{Expe,AccumTest*2} = Recall_O(3);
    
    Resultados_RecallD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_Only_method{Expe,AccumTest*2} = Recall_M(3);
    %%%% FSCORE DETECTION PER CASE    
    Resultados_FscoreD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_SumadoImg{Expe,AccumTest*2} = Fscore_S(3);
    
    Resultados_FscoreD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_OriginalImg{Expe,AccumTest*2} = Fscore_O(3);
    
    Resultados_FscoreD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_Only_method{Expe,AccumTest*2} = Fscore_M(3);
     
    
    
    save([FolderResults,'Results_',char(SelClass),'_Space.mat'],...,
        'Resultados_SumadoImg','Resultados_acc_SVM','Resultados_Only_method','Resultados_OriginalImg',...,
        'Resultados_PrecisionD_SumadoImg',...,
        'Resultados_PrecisionD_OriginalImg',...,
    'Resultados_PrecisionD_Only_method',...,
    'Resultados_RecallD_SumadoImg',...,
    'Resultados_RecallD_OriginalImg',...,
    'Resultados_RecallD_Only_method',...,
    'Resultados_FscoreD_SumadoImg',...,
    'Resultados_FscoreD_OriginalImg',...,
    'Resultados_FscoreD_Only_method')
    
if Expe==46
    ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));    
    imwrite(ImBorde,[FolderImgsSumado,ImgTestS,'.png' ]);
   
   
   
%Imagen borde nucleos detectados en la original
    ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
    imwrite(ImBorde,[FolderImgsOriginal,ImgTestS,'.png' ]);


    
    %%%%% Imagen sin ruido
  imwrite(ImEnhance,[FolderImgsSinRuido,ImgTestS,'.png' ]);
end
 
    
  histo_img = [histo_img;histograms];
  labels_imgTest = [labels_imgTest;labelsAll];

       
end


    end
    end
    AccumTest=0; %variable is reset for the next tissue if is needed
end






