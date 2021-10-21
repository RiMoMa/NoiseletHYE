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
FolderResults = '//mnt/md0/ricardo/NoiseletProject/Results_SVM_2021_B/';
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

%%% Dataset Dir %%%d
FolderIMG = '/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Tissue_images/';
FolderXML = '/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Annotations/';
FolderMasksTest = [FolderIMG,'Masks/'];
FolderTraining = FolderIMG;
FolderTesting = FolderIMG; %cahnge for future datasets
load('ListImgsMonuseg.mat')

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%OPTIONS
      load([FolderResults,'MatrizExperimentosSVM.mat'])%% archivo que contiene los parametros experimentales
      
      
for ClassExp = 2 %realizar experimento por clase
SelClass = clasesMonu(ClassExp)    ;
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

    %%%%%%%%%%%%%
    
    
    %%%%% LEAVE ONE OUT VALIDATION SCHEME %%%%
    
    for OneOut = 1:length(imageForExperiment)
       TestImg =  imageForExperiment(OneOut);
       idxTraining = not(imageForExperiment==TestImg);
       ImgsTraining = imageForExperiment(idxTraining);
       
       %%%%%%% Training %%%%%%%%%
       ListIMGALL = ImgsTraining; %just a change of variable
       ListImgsTrain = ListIMGALL;



    for Expe = 9:16%6
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
%%% VOCABULARY FEATURE EXTRACTION
for Lo = 1:length(ListIMGALL)
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ListIMGALL(Lo)) 
    ImTest = imread([FolderTraining,char(ListIMGALL(Lo)),'.tif']);
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
    ImgTrainName = char(ListImgsTrain(Lo));
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ImgTrainName) 
    
    ImTrain = imread([FolderTraining,ImgTrainName,'.tif']);
    fprintf('Mejorando Imagen\n') 
  %% abrir mascara manual para generar etiquetas
  %  ImGroundT = imread([FolderTraining,ListImgsTrain(Lo).name(1:end-8),'labeled_mask_corrected.png'])>0;
    [binary_maskTest,color_maskTest] = xlmToMask(ImgTrainName,FolderXML,FolderIMG);%% Para
  %  dataset Monuseg
     ImGroundT = binary_maskTest>0;
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
    fprintf('AbriendoImagen\n') 
    ImgTestS = char(TestImg(Lo));
    fprintf('%s\n',ImgTestS) 
    ImTest = imread([FolderTesting,ImgTestS,'.tif']);
    fprintf('Mejorando Imagen\n') 
    MaskName = [FolderMasksTest,ImgTestS,'.png'];
        %% abrir mascara manual para generar etiquetas
tic
    if ~exist(MaskName)
        fprintf('generating image \n')
    [ImGroundT,color_mask]=xlmToMask(ImgTestS,FolderXML,FolderTesting);
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
%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Sin Sumar\n') 
   DiceCoeff = dice(double(MaskEvaluate>0), double(ImManualMask>0)); 
   Resultados_Only_method{Expe,2*OneOut-1} = ListIMGALL(Lo);    
   Resultados_Only_method{Expe,OneOut*2} = DiceCoeff;


   %%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Watershed Original\n') 
   DiceCoeff = dice(double(MaskOriginal>0), double(ImManualMask>0)); 
   Resultados_OriginalImg{Expe,2*OneOut-1} = ImgTestS;    
   Resultados_OriginalImg{Expe,OneOut*2} = DiceCoeff;


   %save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    

   
 %  ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
  % imwrite(ImBorde,[FolderImgsOriginal,ListIMGALL(Lo).name(1:end-4),'.png' ]);

   
   
   
   fprintf('Calculando Metricas Suma Metodos\n') 

   
    DiceCoeff = dice(double(sumaBinary>0), double(ImManualMask>0)); 
    Resultados_SumadoImg{Expe,2*OneOut-1} =ImgTestS;
    Resultados_SumadoImg{Expe,OneOut*2} = DiceCoeff;
    
    
  %  save([FolderResults,'Results_',ResultsName,'_SVM.mat'],'accuaracy','cm')    


    Resultados_acc_SVM{Expe,2*OneOut-1} = ImgTestS;
    Resultados_acc_SVM{Expe,OneOut*2} = accuaracy;

    
    
    
    save([FolderResults,'Results_',char(SelClass),'_Space.mat'],...,
        'Resultados_SumadoImg','Resultados_acc_SVM','Resultados_Only_method','Resultados_OriginalImg')

    %ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));    
    %imwrite(ImBorde,[FolderImgsSumado,ListIMGALL(Lo).name(1:end-4),'.png' ]);
   
   
   
    %%Imagen borde nucleos detectados en la original
 %   ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
 %   imwrite(ImBorde,[FolderImgsOriginal,ImgTestS,'.png' ]);


    
    %%%%% Imagen sin ruido
%  imwrite(ImEnhance,[FolderImgsSinRuido,ImgTestS,'.png' ]);
   
 
    
  histo_img = [histo_img;histograms];
  labels_imgTest = [labels_imgTest;labelsAll];

       
end


    end
    end
end






