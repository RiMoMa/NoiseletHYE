
%pe = pyenv('Version','~/anaconda3/envs/IntelEnv/bin/python',"ExecutionMode","OutOfProcess");
addpath('/home/ricardo/Documents/Doctorado/NucleosMTC/codes/staining_normalization/')
FolderDataset = '/home/ricardo/Documents/Doctorado/NucleosMTC/DatasetTCIA/manual_segmentation_data/';
load('TCIADetail.mat')
ListOfImgCases = DetailData;

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%Folders Output
%%% 
Main = '/home/ricardo/Doctorado/TCIA_PerTissue/';
mkdir(Main)
DatasetFolder = [Main,'Dataset/'];
mkdir(DatasetFolder)  
DatasetFolderAllTisues = [DatasetFolder,'AllTissues_Nuclei_Eval/'];
mkdir(DatasetFolderAllTisues)
DatasetFolderAllTisues_Masks = [DatasetFolderAllTisues,'/Masks/'];
mkdir(DatasetFolderAllTisues_Masks)

DatasetFolderAllTisues_Imgs = [DatasetFolderAllTisues,'/Imgs/'];
mkdir(DatasetFolderAllTisues_Imgs)


DatasetFolderAllTisues_Imgs_Seg = [DatasetFolderAllTisues,'/Imgs_Seg/'];
mkdir(DatasetFolderAllTisues_Imgs_Seg)

DatasetFolderAllTisues_Imgs_SegDual = [DatasetFolderAllTisues,'/Imgs_Seg_dual/'];
mkdir(DatasetFolderAllTisues_Imgs_SegDual)
      
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
    DatasetFolderAllTisues_Imgs_tissue = [DatasetFolderAllTisues_Imgs,num2str(SelClass),'/'];
    mkdir(DatasetFolderAllTisues_Imgs_tissue)
    DatasetFolderAllTisues_Masks_tissue = [DatasetFolderAllTisues_Masks,num2str(SelClass),'/'];
    mkdir(DatasetFolderAllTisues_Masks_tissue)
    DatasetFolderAllTisues_Imgs_tissue_Seg = [DatasetFolderAllTisues_Imgs_Seg,num2str(SelClass),'/'];
    mkdir(DatasetFolderAllTisues_Imgs_tissue_Seg)
    DatasetFolderAllTisues_Imgs_tissue_SegDual = [DatasetFolderAllTisues_Imgs_SegDual,num2str(SelClass),'/'];
    mkdir(DatasetFolderAllTisues_Imgs_tissue_SegDual)

for n =1:length(allImg)
    fprintf('%d de %d\n',n,length(allImg))

    im = imread([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected.png'])>0*255;
    im2 = imread([FolderDataset,'/',num2str(allImg(n)),'_crop.png']);
    ImBorde = imoverlay(im2,boundarymask(im),[155,243,66]/255);
    A = montage({im2,ImBorde});
    A=A.CData;

    %Rutina copiar Imagen
    copyfile([FolderDataset,'/',num2str(allImg(n)),'_crop.png'],DatasetFolderAllTisues_Imgs_tissue)
    imwrite(im,[DatasetFolderAllTisues_Masks_tissue,num2str(allImg(n)),'_crop.png'])

    imwrite(ImBorde,[DatasetFolderAllTisues_Imgs_tissue_Seg,num2str(allImg(n)),'_crop.png'])
    %copyfile([FolderDataset,'/',num2str(allImg(n)),'_labeled_mask_corrected.png'],DatasetFolderAllTisues_Masks)
    imwrite(A,[DatasetFolderAllTisues_Imgs_tissue_SegDual,num2str(allImg(n)),'_crop.png'])

end




end


