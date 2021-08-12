%addpath C:\Users\daedro\Documents\Development\extern\macenko;
%addpath(genpath('C:\Users\daedro\Documents\tmpColor\CPD2'))

a=imread('03.tif');
b=imread('A08_00Ca.tiff');

h=size(b,1);
w=size(b,2);

Io=240;

A= double(reshape(a,[],3));
B= double(reshape(b,[],3));

ODA = double(-log((A+1)/Io));
ODB = double(-log((B+1)/Io));

% scatter3(A(1:100:end,1), A(1:100:end,2), A(1:100:end,3), [], 'b.');
% hold on;
% scatter3(B(1:100:end,1), B(1:100:end,2), B(1:100:end,3), [], 'g.');

%figure;
%scatter3(ODA(1:100:end,1), ODA(1:100:end,2), ODA(1:100:end,3), [], 'b.');
%hold on;
%scatter3(ODB(1:100:end,1), ODB(1:100:end,2), ODB(1:100:end,3), [], 'g.');

% puntosA=pointCloud(ODA(1:500:end,:));
% puntosB=pointCloud(ODB(1:500:end,:));
% 
% tform=pcregrigid(puntosB, puntosA);
% 
% puntosBall=pointCloud(ODB);
% puntosBrectif=pctransform(puntosBall,tform);
% 
% ODBtrans=puntosBrectif.Location;
% 
% Btrans = Io*exp(-ODBtrans);
% Inorm = reshape(Btrans, h, w, 3);
% Inorm = uint8(Inorm);

a=double(a);
b=double(b);

% [x,y]=meshgrid(1:3:255,1:3:255);
% z=griddata(a(:,1),a(:,2),a(:,3),x,y);
% z2=griddata(b(:,1),b(:,2),b(:,3),x,y);
% 
% 
% figure;
% surf(x,y,z, 'FaceColor', 'blue');
% hold on
% surf(x,y,z2+30, 'FaceColor', 'green');
% %surf(x,y,z2+30, 'FaceColor', 'green');
% 


%edgePointsA=[x(gradient(isnan(z))), y(gradient(isnan(z))), z(gradient(isnan(z)))];
%edgePointsB=[x(gradient(isnan(z2))), y(gradient(isnan(z2))), z2(gradient(isnan(z2)))];

%edgePointsA(isnan(edgePointsA(:,3)),:)=[];
%edgePointsB(isnan(edgePointsB(:,3)),:)=[];




opt.method='rigid'; % use rigid registration
opt.viz=1;          % show every iteration
opt.outliers=0;     % do not assume any noise 

opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

opt.max_it=100;     % max number of iterations
opt.tol=1e-8;       % tolerance


% registering Y to X
idxA=round(linspace(1,size(A,1),5000));
idxB=round(linspace(1,size(B,1),5000));

[Transform, Correspondence]=cpd_register(A(idxA,:), B(idxB,:),opt);
[Transform2, Correspondence2]=cpd_register(ODA(idxA,:), ODB(idxB,:),opt);

%[Transform, Correspondence]=cpd_register(edgePointsA, edgePointsB,opt);

figure,cpd_plot_iter(A(idxA,:), B(idxB,:)); title('Before');
figure,cpd_plot_iter(A(idxA,:), Transform.Y);  title('After registering');

figure,cpd_plot_iter(ODA(idxA,:), ODB(idxB,:)); title('Before -OD');
figure,cpd_plot_iter(ODA(idxA,:), Transform2.Y);  title('After registering -OD');

Bnew=cpd_transform(B, Transform);
bnorm = reshape(Bnew, h, w, 3);

Bnew2=cpd_transform(ODB, Transform2);
Bnew2=Io*exp(-Bnew2);
bnorm2 = reshape(Bnew2, h, w, 3);


figure; subplot(1,4,1); 
imshow(uint8(b)); title('Original'); 

subplot(1,4,2); 
imshow(uint8(bnorm)), title('Normalized-pointcloud');

subplot(1,4,3); 
imshow(normalizeStaining(b)), title('Normalized-macenko');

subplot(1,4,4); 
imshow(uint8(bnorm2)), title('Normalized-pointcloud');
