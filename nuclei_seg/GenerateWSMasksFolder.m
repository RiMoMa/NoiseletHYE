function GenerateWSMasksFolder( inputFolder )
%GENERATEWSMASKSFOLDER Generates the watershed masks for existing images in
%the set folder

%addpath /mnt/pan/Data7/gxc206/lymp/code/nuclei_seg/staining_normalization/
%addpath /mnt/pan/Data7/gxc206/lymp/code/nuclei_seg/veta_watershed/

folderContents = dir([inputFolder '/*/*.jpg']);
numFiles=size(folderContents ,1);

mkdir([inputFolder '_masks/']);

for i=numFiles:-1:1
    [~, filename, ~] = fileparts(folderContents(i).name);
    outputFile=[inputFolder '_masks/' filename '_mask.png'];
    
    
    if ~exist(outputFile,'file')
        imgFile=[folderContents(i).folder,'/',folderContents(i).name];
        
        I=imread(imgFile);    
        M = getWatershedMask(I,true,4,10);
        %M = getWatershedMask(I);
        imwrite(M,outputFile);
        else
            fprintf('alredy exist: %s\n',outputFile)
       
    end
end

end

