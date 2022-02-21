
%pe = pyenv('Version','~/anaconda3/envs/IntelEnv/bin/python',"ExecutionMode","OutOfProcess");
addpath('/home/ricardo/Documents/Doctorado/NucleosMTC/codes/staining_normalization/')
FolderDataset = '/home/ricardo/Documents/Doctorado/NucleosMTC/DatasetTCIA/manual_segmentation_data/';
load('TCIADetail.mat')
ListOfImgCases = DetailData;

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%Folders Output
%%% 
Main = '/home/ricardo/Doctorado/UnetPleomorphism/';
mkdir(Main)
DatasetFolder = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/';
mkdir(DatasetFolder)  
DatasetFolderAllTisues = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/AllTissues_Train/';
mkdir(DatasetFolderAllTisues)
DatasetFolderAllTisues_Masks = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/AllTissues_Train/Masks/';
mkdir(DatasetFolderAllTisues_Masks)
DatasetFolderAllTisues_Imgs = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/AllTissues_Train/Imgs/';
mkdir(DatasetFolderAllTisues_Imgs)

DatasetFolderBreast= '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/Breast_Train/';
mkdir(DatasetFolderBreast)
DatasetFolderBreast_Masks = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/Breast_Train/Masks/';
mkdir(DatasetFolderBreast_Masks)
DatasetFolderBreast_Imgs = '/home/ricardo/Doctorado/UnetPleomorphism/Dataset/Breast_Train/Imgs/';
mkdir(DatasetFolderBreast_Imgs)


      
for ClassExp = 1:length(clasesMonu) %realizar experimento por clase
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
    
          

allImg = imageForExperiment;
%inicializa los vectores de las bases

for n =1:length(allImg)
    fprintf('%d de %d\n',n,length(allImg))
    im = imread([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected.png'])>0*255;

   if ~(SelClass == 'brca')
    %Rutina copiar Imagen
    copyfile([FolderDataset,'/',num2str(allImg(n)),'_crop.png'],DatasetFolderAllTisues_Imgs)
    imwrite(im,[DatasetFolderAllTisues_Masks,num2str(allImg(n)),'_crop.png'])
    %copyfile([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected.png'],DatasetFolderAllTisues_Masks)
   else
    copyfile([FolderDataset,'/',num2str(allImg(n)),'_crop.png'],DatasetFolderBreast_Imgs)
    %copyfile([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected.png'],DatasetFolderBreast_Masks)
    imwrite(im,[DatasetFolderBreast_Masks,num2str(allImg(n)),'_crop.png'])
   end

   %im = imread([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected']);%61_labeled_mask_corrected
    %rutina copiar mascara
    
end
% guarda las bases del dataset



end


