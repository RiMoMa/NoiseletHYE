function [X,coords] = NoiseletsFeatures(im,Scales,WinPlsa)
%%% Sacar Mascaras
addpath('PLSA_TEM_EM')
addpath('FNT')
%%% OPTIONS %%%
cluster='kmeans'; % 'plsa'



%K_clusters = 8;
%[Inorm H E] = normalizeStaining(im);
%[H,E] = ColorSepCimalab(im);
%labim = rgb2lab(H);
%ime = labim(:,:,2);
%%ime = ((ime-min(ime(:)))/(max(ime(:))-min(ime(:))))*255;
%ime=uint8(ime);
%ime = im(:,:,1);
% %%%%%%test 
% I = ime
% thresh = multithresh(I,2);
% valuesMax = [thresh max(I(:))];
% [quant8_I_max, index] = imquantize(I,thresh,valuesMax);
% valuesMin = [min(I(:)) thresh];
% quant8_I_min = valuesMin(index);
% imshowpair(quant8_I_min,quant8_I_max,'montage')
% title('Minimum Interval Value           Maximum Interval Value')
% ime = quant8_I_min;




[aa, bb, cc] = size(im);



for win_size = Scales
    
    [Inorm,H,E] = normalizeStaining(im,im,250,0.05);
    
    %[H,E] = ColorSepCimalab(im);
    %labim = rgb2lab(H);
    ime = double( H(:,:,1));
    ime = (ime-min(ime(:)))/(max(ime(:))-min(ime(:)));
    
    
    %%% Extract Noiselets from tiles
    [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(win_size,ime,aa,bb);
    %Rigth shift
    ime_c = zeros(size(ime));
    ime_c(win_size/2:end,:) = ime(win_size/2:end,:);
    [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(win_size,ime_c,aa,bb);
    
    %Down shift
    ime_c = zeros(size(ime));
    ime_c(:,win_size/2:end) = ime(:,win_size/2:end);
    [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(win_size,ime_c,aa,bb);
    % right and down
    
    ime_c = zeros(size(ime));
    ime_c(win_size/2:end,win_size/2:end) = ime(win_size/2:end,win_size/2:end);
    [pila_patches_orig,coords,lienzo,pila_patches_gray] = ExtractTilesAndNoiselets(win_size,ime_c,aa,bb);
    
        
        
    abs_pila = real(pila_patches_orig); %Se saca la Magnitud
    angle_pila = imag(pila_patches_orig); %SE SACA EL ANGULO, se le su

    %vector a clusterizacion
    X=zeros(size(pila_patches_orig,1),size(pila_patches_orig,2)*2);
    X(:,1:2:end-1)=abs_pila;
    X(:,2:2:end) = angle_pila;
    
end


    
end


