
%pe = pyenv('Version','~/anaconda3/envs/IntelEnv/bin/python',"ExecutionMode","OutOfProcess");
addpath('/home/ricardo/Documents/Doctorado/NucleosMTC/codes/staining_normalization/')
FolderDataset = '/home/ricardo/Documents/Doctorado/NucleosMTC/DatasetTCIA/manual_segmentation_data/';
load('TCIADetail.mat')
ListOfImgCases = DetailData;

clasesMonu=unique([ListOfImgCases{:,2}]);%sacar las clases del dataset

%%%%%%%OPTIONS
     
      
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
    
          



%Normalizar dataset Mitos Atypia

 %ruta de las imagenes
%allImg = dir('/home/ricardo/Documents/Doctorado/NucleosMTC/DatasetTCIA/manual_segmentation_data/*p.png');
allImg = imageForExperiment;
%inicializa los vectores de las bases

BasesHE =zeros(size(allImg,1),6);
BasesMaxC =zeros(size(allImg,1),2);
%%
for n =1:length(allImg)
    fprintf('%d de %d\n',n,length(allImg))
    %abre imagen
   im = imread([FolderDataset,'/',num2str(allImg(n)),'_crop.png']);
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
save(['BasesDatasetTIA_',char(SelClass),'.mat'],'BasesHE','BasesMaxC')

%%% AQUI Iria la clusterización y la seleccion de varias bases

%[cluster_n, relatives, ~,~,idx ] = kmedoids(BasesHE,10,'distance','cosine');
pyrun("from sklearn.mixture import BayesianGaussianMixture")
pyrun("dpgmm = BayesianGaussianMixture(n_components=30,covariance_type='full',weight_concentration_prior_type='dirichlet_process',init_params='random',max_iter=200,random_state=2)")
Dic_model = pyrun("q=dpgmm.fit(X)","q",X=BasesHE);
cluster_n = int64(Dic_model.predict(BasesHE))+1;

BasesHE_cluster = zeros(length(unique(cluster_n)),size(BasesHE,2));
BasesMaxC_cluster = zeros(length(unique(cluster_n)),size(BasesMaxC,2));

for cln = unique(cluster_n)
    lo_cluster = cluster_n==cln;
    bases_cluster_HE = median(BasesHE(lo_cluster,:),1);
    bases_cluster_HE_MAX = median(BasesMaxC(lo_cluster,:),1);

   BasesHE_cluster(cln,:) = bases_cluster_HE;
   BasesMaxC_cluster(cln,:) = bases_cluster_HE_MAX;
end

%idx = pyrunfile("/home/ricardo/PycharmProjects/Normalizacion/ClusterData.py","ClusterAssigned",X=BasesHE)

%MaxC_clusters = BasesMaxC(idx,:);




%calcula la media y genera una base HE haciendo reshape
HENORM = median(BasesHE);
HENORM =reshape(HENORM,3,2);
HEMaxC = median(BasesMaxC)';




%ruta de salida de las imagenes normalizadas
FolderData = '/home/ricardo/Documents/Doctorado/DatasetTCIA_NORM_new/';
mkdir(FolderData)
FolderComposed = ['/home/ricardo/Documents/Doctorado/DatasetTCIA_Norm_groups/',char(SelClass),'/'];
mkdir(FolderComposed)

%normaliza todo el dataset Con base a la media seleccionada HENORM
parfor n =1:length(allImg)
   fprintf('%d de %d\n',n,length(allImg))
   im = imread([FolderDataset,'/',num2str(allImg(n)),'_crop.png']);
   %%% TAKE the values of the representant to normalize
   takefrom = cluster_n(n);
   mkdir([FolderComposed,num2str(takefrom),'/'])
   HENORM = BasesHE_cluster(takefrom,:);
   HENORM =reshape(HENORM,3,2);
   HEMaxC = BasesMaxC_cluster(takefrom,:)';
   [Inorm H E HE] = normalizeStaining(im,im, 250, 0.05, 1, HENORM,HEMaxC);
   %clase  = allImg(n).folder;
   Compose2 = montage({im,Inorm});
   Compose2 = Compose2.CData;
   name = num2str(allImg(n));%.name;
 %  clase=clase(end-6:end);
   name = [name,'_crop.png'];
   %Ruta donde guarda la imagen normalizada nueva
   rutaIMG = [FolderData,'/',name];
   imwrite(Inorm,rutaIMG);   
   
   rutaCompo = [FolderComposed,'/',num2str(takefrom),'/',name];
   imwrite(Compose2,rutaCompo);
    
end

end