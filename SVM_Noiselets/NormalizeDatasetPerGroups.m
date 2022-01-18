addpath('/home/ricardo/Documents/Doctorado/NucleosMTC/codes/staining_normalization/')

%Normalizar dataset Mitos Atypia

 %ruta de las imagenes
allImg = dir('/home/ricardo/Documents/Doctorado/NucleosMTC/DatasetTCIA/manual_segmentation_data/*p.png');
%inicializa los vectores de las bases
BasesHE =zeros(size(allImg,1),6);
BasesMaxC =zeros(size(allImg,1),2);
%%
for n =1:length(allImg)
    fprintf('%d de %d\n',n,length(allImg))
    %abre imagen
   im = imread([allImg(n).folder,'/',allImg(n).name]);
   %calcula las bases usando masenko
   %(I,ImO, Io, beta, alpha, HERef, maxCRef)
   [Inorm H E HE maxC] = normalizeStaining(im,im,250,0.05);
   HE=HE(:)';
   maxC = maxC(:)';
   %Guarda las bases encontradas en el vector
   BasesHE(n,:)=HE;
   BasesMaxC(n,:)=maxC;
    
end
% guarda las bases del dataset
save('BasesDatasetMitos.mat','BasesHE','BasesMaxC')

%%% AQUI Iria la clusterización y la seleccion de varias bases

[cluster_n, relatives, ~,~,idx ] = kmedoids(BasesHE,10,'distance','cosine');

MaxC_clusters = BasesMaxC(idx,:);




%calcula la media y genera una base HE haciendo reshape
HENORM = median(BasesHE);
HENORM =reshape(HENORM,3,2);
HEMaxC = median(BasesMaxC)';




%ruta de salida de las imagenes normalizadas
FolderData = '/home/ricardo/Documents/Doctorado/DatasetTCIA_NORM/';
mkdir(FolderData)
FolderComposed = '/home/ricardo/Documents/Doctorado/DatasetTCIA_Norm_groups/';
mkdir(FolderComposed)

%normaliza todo el dataset Con base a la media seleccionada HENORM
parfor n =1:length(allImg)
   fprintf('%d de %d\n',n,length(allImg))
   im = imread([allImg(n).folder,'/',allImg(n).name]);
   %%% TAKE the values of the representant to normalize
   takefrom = cluster_n(n);
   mkdir([FolderComposed,num2str(takefrom),'/'])
   HENORM = BasesHE(idx(takefrom),:);
   HENORM =reshape(HENORM,3,2);
   HEMaxC = BasesMaxC(idx(takefrom),:)';
   [Inorm H E HE] = normalizeStaining(im,im, 250, 0.05, 1, HENORM,HEMaxC);
   %clase  = allImg(n).folder;
   Compose2 = montage({im,Inorm});
   Compose2 = Compose2.CData;
   name = allImg(n).name;
 %  clase=clase(end-6:end);
   name = [name(1:end-4),'.png'];
   %Ruta donde guarda la imagen normalizada nueva
   rutaIMG = [FolderData,'/',name];
   imwrite(Inorm,rutaIMG);   
   
   rutaCompo = [FolderComposed,'/',num2str(takefrom),'/',name];
   imwrite(Compose2,rutaCompo);
    
end