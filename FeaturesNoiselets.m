function [Xnucleos,Xfondo] = FeaturesNoiselets (ImgH,ImgGroundT,nSample)
%nSample = 4;
Window_size = 2^nSample;
imgY = size(ImgH,1);
imgX = size(ImgH,2);
porcionX = imgX/Window_size;
porcionY = imgY/Window_size;

ImgH=rgb2gray(ImgH);
%matriz Noiselets
NoiseletMatrix = generate_noiselet_recur(Window_size);
%imagenFondo
imFondo = (im2double(ImgH).*not(double(ImgGroundT>0)));
imFondo = imFondo+double(ImgGroundT>0);
%imagenNucleos
imNucleos = (im2double(ImgH).*double(ImgGroundT>0));

%cuadricular
Xnucleos = [];
Xfondo =[];
for X =1:porcionX
for Y=1:porcionY
    imgFondoAux = imFondo((Y-1)*Window_size+1:(Y)*Window_size,(X-1)*Window_size+1:(X)*Window_size);
    imgNucleosAux = imNucleos((Y-1)*Window_size+1:(Y)*Window_size,(X-1)*Window_size+1:(X)*Window_size);
    
   if ~all(all(imgNucleosAux==0)==1)
    imgNucleosAuxNSL =  imgNucleosAux.*NoiseletMatrix;
    Xnucleos = [Xnucleos;imgNucleosAuxNSL(:)'];
   end
    imgFondoAuxNSL =  imgFondoAux.*NoiseletMatrix;
    Xfondo = [Xfondo;imgFondoAuxNSL(:)'];

    
end
end
%aplicar noiselets
%etiquetar
%aplicar noiselets
%etiquetar
end