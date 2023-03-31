%%% REorganizing Experiments pipeline Noiselets 

clear all

UsedDistance = 'cityblock';
%%%% Libraries %%%%%%%%%%%
addpath('/mnt/md0/ricardo/CodesInternship/GITP/')
%addpath('/mnt/storage/ramoncayomar/Codes/CodesInternship/GITP/')
%addpath('/home/ricardo/Documents/Doctorado/SYNC/all_codes/CodesInternship/GITP/')

addpath('nuclei_seg/staining_normalization/')
addpath('SVM_Noiselets/')
addpath('FNT')
addpath('PLSA_TEM_EM/')
addpath('nuclei_seg')
addpath('nuclei_seg/veta_watershed/')
%%%%%%%%%%%%%%%%%

%%%%%%%%% Folder Results %%%%%%%%
FolderResults = '/mnt/md0/ricardo/NoiseletProject/ResultsNoiselets_Original_Complete_noOverlap_Monuseg_TestSeminar/'; %folder Out'
%FolderResults =     '/home/ricardo/Documents/Doctorado/ResultsNoiselets_Original_Complete_noOverlap_Monuseg/'; %folder Out
%FolderResults =     '/mnt/storage/ramoncayomar/ResultsNoiselets_Original_Complete_noOverlap_Monuseg/'; %folder Out
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
FolderCoefficients  = [FolderResults,'AllCoefficients/']; %% borderline groundtruth
mkdir(FolderCoefficients);
FolderHistogramsTrain  = [FolderResults,'AllTrainHists/']; %% borderline groundtruth
mkdir(FolderHistogramsTrain);
FolderHistogramsTest  = [FolderResults,'AllTestHists/']; %% borderline groundtruth
mkdir(FolderHistogramsTest);
%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Dataset Dir %%%
%FolderIMG = '/mnt/md0/Histopathology/Datasets/TCIA/';
%FolderIMG = '/mnt/storage/ramoncayomar/Datasets/TCIA/';   
%FolderIMG = '/home/ricardo/Documents/Doctorado/DatasetTCIA_NORM/';
%FolderXML =
%'/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Annotations/'; %this is for monuseg dataset
FolderIMG = '/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Tissue_images/';
FolderXML = '/mnt/md0/Histopathology/Datasets/DatasetMonusegC/Annotations/';


FolderMasksTest = [FolderIMG,'Masks/'];
FolderTraining = FolderIMG;
FolderTesting = FolderIMG; %change for future datasets
load('ListImgsMonuseg.mat')
%ListOfImgCases = DetailData; %TCIA ONLY

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%OPTIONS
load([FolderResults,'MatrizExperimentosSVM_2022.mat'])%% archivo que contiene los parametros experimentales
      

for ClassExp = 3%:length(clasesMonu) %realizar experimento por clase
    imageForExperiment={};
    SelClass = clasesMonu(ClassExp) ;
    fprintf('Evaluating for %s tissue\n',SelClass)
    Variable_Data_arragment = strcat(FolderResults,'Data_Array_',SelClass,'.mat');
   if ~exist(Variable_Data_arragment,'file')  
    count=0; 
    for rn = 1:length(ListOfImgCases)%take out list per tissue selection
       ItClass = ListOfImgCases{rn,2}==SelClass; % logical existence of a determine class
       
       if ItClass
          count=count+1; 
          imageForExperiment{count} = ListOfImgCases{rn,1} ;
       end
    end
    
    
     
    %%%%%%%%%%%%% K-Fold Partition
%  KFold = 10;
    idxForExperiment = 1:length(imageForExperiment);
    cv = cvpartition(idxForExperiment,'leaveout');
    %cv = cvpartition(imageForExperiment,'Holdout',0.3)
    save(Variable_Data_arragment,'cv','idxForExperiment','imageForExperiment')
    else 
        
        load(Variable_Data_arragment)
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
    
   

   
   
for Expe = 1%:length(ExperimentosNoiselets)
        % Load parameters for each experiment
        Scales = 16%[ExperimentosNoiselets{Expe,2}];
        WinPlsa = 16%ExperimentosNoiselets{Expe,3};
        ResultsName = ExperimentosNoiselets{Expe,1}{1};
        K_clusters = ExperimentosNoiselets{Expe,4};
        fprintf('Testing: %s \n',ResultsName)
        ExpeName = strcat('_case_',SelClass,'_Exp_',num2str(Expe));
        ExpeName2 = strcat('_case_',SelClass,'_Scale_',num2str(Scales),'_WinPlsa_',num2str(WinPlsa));
        ExpeName2b = strcat('_case_',SelClass,'_Scale_',num2str(Scales));
        ExpeName3 = strcat('_case_',SelClass,'_Scale_',num2str(Scales),'_WinPlsa_',num2str(WinPlsa),'_clusters_',num2str(K_clusters));
 
        AccumTest = 0;     
 for OneOut = 1:cv.NumTestSets
       TestImg = imageForExperiment(cv.test(OneOut));
       idxTraining = cv.training(OneOut);
       ImgsTraining = imageForExperiment(idxTraining);
       
       %%%%%%% Training %%%%%%%%%
       ListIMGALLT = ImgsTraining; %just a change of variable
       ListImgsTrain = ListIMGALLT;
       
     Variable_Vocab_cases = strcat(FolderResults,'VocabList_',ExpeName2,'_Fold_',num2str(OneOut),'.mat');

       if ~exist(Variable_Vocab_cases,'file')
        %VocabCases = randperm(length(ListIMGALLT),round(length(ListIMGALLT)*0.4));
        VocabCases = randperm(length(ListIMGALLT),length(ListIMGALLT));
        save(Variable_Vocab_cases,'VocabCases')
       else
           load(Variable_Vocab_cases)
       end
       
	%% Feature Extraction For dictionary Building
    allFeatures = [];
    allLabels = [];
    ListIMGALL = ListIMGALLT(VocabCases);

    for Lo = 1:length(ListIMGALL)
    ImageName = num2str(ListIMGALL{Lo});
    %idx = findstr(ImageName,'.png')
    fprintf('AbriendoImagen\n') 
    fprintf('%s\n',ImageName)
    PathCoefficients = strcat(FolderCoefficients ,ImageName, ExpeName2b,'.mat');
    MaskName = [FolderMasksTest,ImageName,'.mat'];
      if ~exist(MaskName,'file')
  %% abrir mascara manual para generar etiquetas
  % ImGroundT = imread([FolderTraining,ImgTrainName,'_labeled_mask_corrected.png'])>0;
   [ImGroundTE,color_maskTest] = xlmToMask(ImgTrainName,FolderXML,FolderIMG);%% Para
      save(MaskName,'ImGroundTE')
    
    else
       ImGroundTE = load(MaskName).ImGroundTE;       
        
    end
    if ~exist(PathCoefficients)
    ImTest = imread([FolderTraining,ImageName,'.tif']);
    fprintf('Features for vocab building\n') 
    %ToDo: COncatenate img features
    [aa,bb,cc] = size(ImTest);
    GroundT = ImGroundTE>0;
    [Inorm,H,E] = normalizeStaining(ImTest,ImTest,220,0.06);
        ime = double(H(:,:,1));
    ime=(ime-min(ime(:)))/(max(ime(:))-min(ime(:)));
    [pila_patches_orig,coords,lienzo,pila_patches_gray,labels] = ExtractTilesAndNoiseletsLabels_scales(Scales,ime,aa,bb,GroundT);
   % [XFeatures,CoordIMG] = NoiseletsFeatures(ImTest,Scales);
    save(PathCoefficients,'pila_patches_orig','labels')
    else
        load(PathCoefficients)
    end
    allFeatures = [allFeatures;pila_patches_orig];
    allLabels = [allLabels;labels];    
    end  
    
% idxF = randperm(10000)    ;
% allFeatures = allFeatures(idxF,:);
% allLabels=allLabels(idxF);
% abs_pila = real(allFeatures); %Se saca la Magnitud
% angle_pila = imag(allFeatures); %SE SACA EL ANGULO, se le su
% %vector a clusterizacion
% X=zeros(size(allFeatures,1),size(allFeatures,2)*2);
% X(:,1:2:end-1)=abs_pila;
% X(:,2:2:end) = angle_pila;
% XFeatures = X;
% Y1 =tsne(XFeatures,'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('All Features all vector')
% Y1 =tsne(abs_pila,'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('all features Magnitude')
% Y1 =tsne(angle_pila,'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('all F imag')
% Y1 =tsne(abs_pila(:,1:256),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('1scale real')
% Y1 =tsne(angle_pila(:,1:256),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('1scale imag')
% Y1 =tsne(abs_pila(:,257:512),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('2scale real')
% Y1 =tsne(angle_pila(:,257:512),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('2scale imag')
% Y1 =tsne(abs_pila(:,513:256*3),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('3scale real')
% Y1 =tsne(angle_pila(:,513:256*3),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('3 scale imag')
% Y1 =tsne(abs_pila(:,256*3+1:1024),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('4 scale real')
% Y1 =tsne(angle_pila(:,256*3+1:1024),'Distance','cityblock');
% figure;gscatter(Y1(:,1),Y1(:,2),allLabels);title('4scale imag')
%% build the dictionary using kmeans

NameVocab = strcat(FolderCoefficients,'Vocab_',ExpeName,'_Fold_',num2str(OneOut),'.mat');
if ~exist(NameVocab,'file')
    
a=allFeatures(:,1:256);
b=allFeatures(:,257:512);
c=allFeatures(:,513:256*3);
d=allFeatures(:,256*3+1:1024);

edges = [-2:4/100:2];
allF=zeros(size(allFeatures,1),(size(edges,2)-1)*8);

parfor n=1:size(allFeatures,1)
F1=histogram(real(a(n,:)),'BinEdges',edges);
F1=F1.Values;
F2=histogram(real(b(n,:)),'BinEdges',edges);
F2=F2.Values;
F3=histogram(real(c(n,:)),'BinEdges',edges);
F3=F3.Values;
F4=histogram(real(d(n,:)),'BinEdges',edges);
F4=F4.Values;
F1b=histogram(imag(a(n,:)),'BinEdges',edges);
F1b=F1b.Values;
F2b=histogram(imag(b(n,:)),'BinEdges',edges);
F2b=F2b.Values;
F3b=histogram(imag(c(n,:)),'BinEdges',edges);
F3b=F3b.Values;
F4b=histogram(imag(d(n,:)),'BinEdges',edges);
F4b=F4b.Values;
allF(n,:) =[F1,F1b,F2,F2b,F3,F3b,F4,F4b];
end


% 
     save(NameVocab,'allF','allLabels')
     else
         load(NameVocab)
     end


%     %Build histograms and extract labels
%     histo_img=[];
%     labels_img=[];
% for Lo = 1:length(ListImgsTrain)
%     ImgTrainName = num2str(ListImgsTrain{Lo});
%     fprintf('AbriendoImagen\n') 
%     fprintf('%s\n',ImgTrainName) 
%     
%     ImTrain = imread([FolderTraining,ImgTrainName,'.tif']);%Todo: Check
%     fprintf('Mejorando Imagen\n') 
%    MaskName = [FolderMasksTest,ImgTrainName,'.mat'];
%   if ~exist(MaskName,'file')
%   %% abrir mascara manual para generar etiquetas
%   % ImGroundT = imread([FolderTraining,ImgTrainName,'_labeled_mask_corrected.png'])>0;
%    [ImGroundTE,color_maskTest] = xlmToMask(ImgTrainName,FolderXML,FolderIMG);%% Para
%       save(MaskName,'ImGroundTE')
%     
%     else
%        ImGroundTE = load(MaskName).ImGroundTE;       
%         
%     end
%   %dataset Monuseg
%      ImGroundT = ImGroundTE>0;
%      PathTrainHistImg = strcat(FolderHistogramsTrain,ImgTrainName,'HistoExp_',ExpeName,'_Fold_',num2str(OneOut),'.mat');
%      if ~exist(PathTrainHistImg,'file')
%       [labelsAll,histograms] = NoiseletsPLSAHistogramImg(ImTrain,Scales,WinPlsa,ImGroundT,vocab,[],UsedDistance);
%       save(PathTrainHistImg,'labelsAll','histograms')
%      else
%          load(PathTrainHistImg)
%      end
%        histo_img = [histo_img;histograms];
%        labels_img = [labels_img;labelsAll];
%        
% end



%% train a svm
X=allF;
Y=allLabels;


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
    ImgTestS = num2str(TestImg{Lo});
    fprintf('%s\n',ImgTestS) 
    ImTest = imread([FolderTesting,ImgTestS,'.tif']);
    fprintf('Mejorando Imagen\n') 
    MaskName = [FolderMasksTest,ImgTestS,'.mat'];
        %% abrir mascara manual para generar etiquetas
tic
    if ~exist(MaskName,'file')
        fprintf('generating image \n')
    %ImGroundTE = imread([FolderTraining,ImgTestS,'_labeled_mask_corrected.png']);
  %  ImGroundT = ImGroundTE >0; 
        [ImGroundTE,color_mask]=xlmToMask(ImgTestS,FolderXML,FolderTesting);
    
    save(MaskName,'ImGroundTE')
    
    else
       ImGroundTE = load(MaskName).ImGroundTE;       
        
    end
    ImGroundT = ImGroundTE >0;
    ImManualMask=ImGroundT;
toc
  % [binary_maskTest,color_maskTest] =
  % xlmToMask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG);%% Para
  % dataset Monuseg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% Removing Noise of the classified histograms%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         
       try
              
          PathTestHistImg = strcat(FolderHistogramsTest,ImgTestS,'HistoExp_',ExpeName,'_Fold_',num2str(OneOut),'.mat');
        if ~exist(PathTestHistImg,'file')
            [histograms,ImEnhance,ImNormNoiselet] = NoiseletsPLSAHistogramImg_scales(ImTest,Scales,WinPlsa,ImGroundT,ClassModel);
            save(PathTestHistImg,'histograms','ImEnhance','ImNormNoiselet')
         else
         load(PathTestHistImg)
            end
  
         %%%% SACAR RESULTADOS GENERALES DE LA CLASSIFICACION

        %% Evaluate SVM
        Y_predict=predict(ClassTreeEns,histograms);

        %Results



        %[c,cm,ind,per] = confusion(labelsAll',Y_predict');
        %accuaracy = sum(diag(cm))/sum(cm(:));
        %fprintf('accuaracy:%d\n',accuaracy)

    
    catch
        ImEnhance = ImTest;
        accuaracy = 0; % this value was not computed
        fprintf('No se pudo encontrar señal\n')
    end


%%%%Sacando Mascaras
   fprintf('Sacando Mascara\n') 
   MaskEvaluate = getWatershedMask(ImEnhance);
   MaskOriginal = getWatershedMask(ImTest);
   fprintf('Suma de las mascaras \n')
   sumaBinary = or(MaskOriginal,MaskEvaluate);
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% METRICS %%%%%%%%%%%
   fprintf('Calculando Metricas Sin Sumar\n') 
   DiceCoeff = dice(double(MaskEvaluate>0), double(ImManualMask>0)); 
   % DETECTION METRIC
   Resultados_Only_method{Expe,2*AccumTest-1} = ImgTestS;    
   Resultados_Only_method{Expe,AccumTest*2} = DiceCoeff;
   
%%%%%%%%%%%%%%%%%%
   fprintf('Calculando Metricas Watershed Original\n') 
   DiceCoeff = dice(double(MaskOriginal>0), double(ImManualMask>0)); 
   Resultados_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;    
   Resultados_OriginalImg{Expe,AccumTest*2} = DiceCoeff;
    
   
   fprintf('Calculando Metricas Suma Metodos\n') 
    DiceCoeff = dice(double(sumaBinary>0), double(ImManualMask>0)); 
    Resultados_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_SumadoImg{Expe,AccumTest*2} = DiceCoeff;
 
    %Resultados_acc_SVM{Expe,2*AccumTest-1} = ImgTestS;
    %Resultados_acc_SVM{Expe,AccumTest*2} = accuaracy;
    
       %%% Metrics for nuclei Detection
 
[Precision_S,Recall_S,Fscore_S] = Nuclei_centroidsBased_metric(ImGroundTE, sumaBinary);
[Precision_O,Recall_O,Fscore_O] = Nuclei_centroidsBased_metric(ImGroundTE, MaskOriginal);
[Precision_M,Recall_M,Fscore_M] = Nuclei_centroidsBased_metric(ImGroundTE, MaskEvaluate);

       %%%% PRECISION DETECTION PER CASE
    Resultados_PrecisionD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_SumadoImg{Expe,AccumTest*2} = Precision_S(4);
    
    Resultados_PrecisionD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_OriginalImg{Expe,AccumTest*2} = Precision_O(4);
    
    Resultados_PrecisionD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_PrecisionD_Only_method{Expe,AccumTest*2} = Precision_M(4);
        
    %%%% RECALL DETECTION PER CASE    
    Resultados_RecallD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_SumadoImg{Expe,AccumTest*2} = Recall_S(4);
    
    Resultados_RecallD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_OriginalImg{Expe,AccumTest*2} = Recall_O(4);
    
    Resultados_RecallD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_RecallD_Only_method{Expe,AccumTest*2} = Recall_M(4);
    %%%% FSCORE DETECTION PER CASE    
    Resultados_FscoreD_SumadoImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_SumadoImg{Expe,AccumTest*2} = Fscore_S(4);
    
    Resultados_FscoreD_OriginalImg{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_OriginalImg{Expe,AccumTest*2} = Fscore_O(4);
    
    Resultados_FscoreD_Only_method{Expe,2*AccumTest-1} = ImgTestS;
    Resultados_FscoreD_Only_method{Expe,AccumTest*2} = Fscore_M(4);
     
    
    
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
    
if Expe==20
    ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));    
    imwrite(ImBorde,[FolderImgsSumado,ImgTestS,'.png' ]);
   
   
   
%    Imagen borde nucleos detectados en la original
    ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
    imwrite(ImBorde,[FolderImgsOriginal,ImgTestS,'.png' ]);


    
    %%%%% Imagen sin ruido
  imwrite(ImEnhance,[FolderImgsSinRuido,ImgTestS,'.png' ]);
end
 
    
%  histo_img = [histo_img;histograms];
%  labels_imgTest = [labels_imgTest;labelsAll];

       
end


    end
    end
    AccumTest=0; %variable is reset for the next tissue if is needed
end

    
 


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 