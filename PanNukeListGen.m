ListRoutes = {'~/Doctorado/PanNuke/PanNuke_Dataset/Fold1/DatasetImgs/*/*.png';
'~/Doctorado/PanNuke/PanNuke_Dataset/Fold2/DatasetImgs/*/*.png';
'~/Doctorado/PanNuke/PanNuke_Dataset/Fold3/DatasetImgs/*/*.png'};

AllNames={};
AllClasses = {};
AllPaths={};
AllTrainStage={};
for z=1:length(ListRoutes)


Fold1=dir(ListRoutes{z});


Names={};
Classes = {};
Paths={};
TrainStage=[];
fprintf('Imagenes %d\n',length(Fold1))
for n=1:length(Fold1)
    File = Fold1(n).name;
    idx = strfind(File,'.png')-1;
    NameIm = File(1:end-idx);
    Folder =Fold1(n).folder;
    idxType = strfind(Folder,'/');
    idxType = idxType(end)+1;
    Type = Folder(idxType:end);
    pathComplete = [Folder,'/',File];

    Paths{n,1}=convertCharsToStrings(pathComplete);
    Classes{n,1}=convertCharsToStrings(Type);
    Names{n,1}=convertCharsToStrings(NameIm);
    TrainStage{n,1}=z;

end

AllNames=[AllNames;Names];
AllClasses = [AllClasses;Classes];
AllPaths=[AllPaths;Paths];
AllTrainStage=[AllTrainStage;TrainStage];

end
ListOfImgCases = [AllNames,AllClasses,AllPaths,AllTrainStage];
save('/home/ricardo/Documents/Doctorado/NucleosMTC/codes/ListImgsPanNuke.mat','ListOfImgCases')

