function [Xf,imAsoc] = FeaturesNoiseletsTest (ImgH,nSample)
%nSample = 4;
Window_size = 2^nSample;
imgY = size(ImgH,1);
imgX = size(ImgH,2);
porcionX = imgX/Window_size;
porcionY = imgY/Window_size;
%matriz Noiselets
NoiseletMatrix = generate_noiselet_recur(Window_size);

ImgH = im2double(rgb2gray(ImgH));

%cuadricular
Xf = [];
imAsoc = [];
for X =1:porcionX
for Y=1:porcionY
    imAux = ImgH((Y-1)*Window_size+1:(Y)*Window_size,(X-1)*Window_size+1:(X)*Window_size);
       
 
    imgAuxNSL =  imAux.*NoiseletMatrix;
    Xf = [Xf;imgAuxNSL(:)'];
    imAsoc = [imAsoc;imAux(:)'];
    
end
end
%aplicar noiselets
%etiquetar
%aplicar noiselets
%etiquetar
end