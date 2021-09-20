%%%%%%%%%%%%
%Using Noiselets to improve Nuclei Segmentation methods
%@version 2.0
%@author M.Sc. Ricardo Moncayo <ramoncayomar@unal.edu.co>
%Ph.D. Eduardo Romero
%%%%%%%%%%%%



addpath('/mnt/md0/ricardo/CodesInternship/GITP/')
addpath('nuclei_seg/staining_normalization/')
addpath('FNT')
addpath('PLSA_TEM_EM/')
addpath('nuclei_seg')
addpath('nuclei_seg/veta_watershed/')

%%% Dataset Dir %%%
FolderIMG = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Tissue_images/';
FolderXML = '/mnt/md0/Histopathology/Datasets/MoNuSeg/Annotations/';


%% Folder Results
FolderResults = '//mnt/md0/ricardo/NoiseletProject/ResultsMonuseg2021/';
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
load([FolderResults,'MatrizExperimentos3.mat'])%% archivo que contiene los parametros experimentales
for Expe = 97:length(ExperimentosNoiselets)
Scales = [ExperimentosNoiselets{Expe,2}];
WinPlsa = ExperimentosNoiselets{Expe,3};
ResultsName = ExperimentosNoiselets{Expe,1}{1};
K_clusters =ExperimentosNoiselets{Expe,4};
fprintf('Testing: %s \n',ResultsName)
%%%%%%%Training
ListIMGALL = dir(strcat(FolderIMG,'*.tif'));
Resultados = cell(length(ListIMGALL),5);
ResultadosSumado =cell(length(ListIMGALL),5);
ResultadosOriginal=cell(length(ListIMGALL),5);
%%leave one out
for Lo = 1:length(ListIMGALL)
fprintf('AbriendoImagen\n') 
fprintf('%s\n',ListIMGALL(Lo).name) 
    ImTest = imread([FolderIMG,ListIMGALL(Lo).name]);
    fprintf('Mejorando Imagen\n') 
    [ImEnhance,~] = MascaraDeNoiseletsPLSA(ImTest,Scales,WinPlsa,K_clusters);
    fprintf('Sacando Mascara\n') 
    MaskEvaluate = getWatershedMask(ImEnhance);
    imwrite(ImEnhance,[FolderImgsSinRuido,ListIMGALL(Lo).name(1:end-4),'.png' ]);

fprintf('Calculando Metricas\n') 
[binary_maskTest,color_maskTest,JAI,DiceCoeff,JAC] = he_to_binary_mask(ListIMGALL(Lo).name(1:end-4),FolderXML,FolderIMG,MaskEvaluate);
DiceCoeff = dice(double(MaskEvaluate>0), double(binary_maskTest>0)); 
JAC = jaccard(double(MaskEvaluate>0), double(binary_maskTest>0)); 
Resultados{Lo,1} = ListIMGALL(Lo).name(1:end-4);
Resultados{Lo,2} = 0;
Resultados{Lo,3} = DiceCoeff;
Resultados{Lo,4} = JAC;
%aji3 = Aggregated_Jaccard_Index_v1_0(binary_maskTest,bwlabel(MaskEvaluate));
aji3=0;
Resultados{Lo,5} = aji3;

save([FolderResults,'Results_',ResultsName,'_Method.mat'],'Resultados')    
ImBorde = imoverlay(ImTest,boundarymask(MaskEvaluate));
imwrite(ImBorde,[FolderImgsMethod,ListIMGALL(Lo).name(1:end-4),'.png' ]);

fprintf('Mascara Metodo Original \n') 
MaskOriginal = getWatershedMask(ImTest);
fprintf('Calculando Metricas\n') 

DiceCoeff = dice(double(MaskOriginal>0), double(binary_maskTest>0)); 
JAC = jaccard(double(MaskOriginal>0), double(binary_maskTest>0)); 

ResultadosOriginal{Lo,1} = ListIMGALL(Lo).name(1:end-4);
ResultadosOriginal{Lo,2} = 0;
ResultadosOriginal{Lo,3} = DiceCoeff;
ResultadosOriginal{Lo,4} = JAC;
%aji1 = Aggregated_Jaccard_Index_v1_0(binary_maskTest,bwlabel(MaskOriginal));
aji1=0;
ResultadosOriginal{Lo,5} = aji1;
save([FolderResults,'Results_',ResultsName,'Original.mat'],'ResultadosOriginal')    
ImBorde = imoverlay(ImTest,boundarymask(MaskOriginal));
imwrite(ImBorde,[FolderImgsOriginal,ListIMGALL(Lo).name(1:end-4),'.png' ]);


fprintf('Suma de las mascaras \n')
sumaBinary = or(MaskOriginal,MaskEvaluate) ;
DiceCoeff = dice(double(sumaBinary>0), double(binary_maskTest>0)); 
JAC = jaccard(double(sumaBinary>0), double(binary_maskTest>0)); 
ResultadosSumado{Lo,1} = ListIMGALL(Lo).name(1:end-4);
ResultadosSumado{Lo,2} = 0;
ResultadosSumado{Lo,3} = DiceCoeff;
ResultadosSumado{Lo,4} = JAC;
%aji2 = Aggregated_Jaccard_Index_v1_0(binary_maskTest,bwlabel(or(MaskOriginal,MaskEvaluate)));
aji2=0;
ResultadosSumado{Lo,5} = aji2;

save([FolderResults,'Results_',ResultsName,'Sumados.mat'],'ResultadosSumado')    
ImBorde = imoverlay(ImTest,boundarymask(or(MaskOriginal,MaskEvaluate)));
imwrite(ImBorde,[FolderImgsSumado,ListIMGALL(Lo).name(1:end-4),'.png' ]);

ImBorde = imoverlay(ImTest,boundarymask(binary_maskTest>0));
imwrite(ImBorde,[FolderImgsGroundTruth,ListIMGALL(Lo).name(1:end-4),'.png' ])
fprintf('Resultados Para:\n%s\n',ListIMGALL(Lo).name)
fprintf('Original:%d\n',aji1)
fprintf('Metodo:%d\n',aji3)
fprintf('Suma:%d\n',aji2)
fprintf('Diferencia:%d\n',aji1-aji3)

end

end

%sacar cuales mejoraron y cuales empeoraron
FolderImgsMejoran  = [FolderResults,'ImgsMejoran/'];
mkdir(FolderImgsMejoran);
FolderImgsEmpeoran = [FolderResults,'ImgsEmpeoran/'];
mkdir(FolderImgsEmpeoran);
Mejoraron = {};
empeoraron = {};

CompararAJI = [[ResultadosOriginal{:,3}];[ResultadosSumado{:,3}]];
indexMejoraron = CompararAJI(1,:)<CompararAJI(2,:);

for id = 1:length(ind exMejoraron)
    imName = [Resultados{id,1},'.png'];
    ImPathOriginalSegmentation = [FolderImgsOriginal,imName];
    ImPathMethodSegmentation = [FolderImgsMethod,imName];
    ImPathGT = [FolderImgsGroundTruth,imName];
    ImPathSumado = [FolderImgsSumado,imName];
    A = montage({ImPathOriginalSegmentation,ImPathMethodSegmentation,ImPathSumado,ImPathGT},'BackgroundColor','black','BorderSize',[5,5],'ThumbnailSize',[] );
    A = A.CData;
    close all
    if indexMejoraron(id) == 1
        imwrite(A,[FolderImgsMejoran,imName]);
    else
        imwrite(A,[FolderImgsEmpeoran,imName]);
    end
        
end

CompararAJI={};
CompararAJI(:,1) = ResultadosSumado(:,1);
CompararAJI(:,2) = ResultadosOriginal(:,3);
CompararAJI(:,3) = num2cell([double(indexMejoraron)]);
CompararAJI(:,4) = ResultadosSumado(:,3);
CompararAJI(:,5) = num2cell([ResultadosOriginal{:,3}]-[ResultadosSumado{:,3}]);%diferencias

[d,idMaxDiff] = max([CompararAJI{:,3}]);




