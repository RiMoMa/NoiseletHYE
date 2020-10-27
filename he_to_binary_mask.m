% A function to read in H&E image and xml file containing annotations
% Gives the binary and colored maps based on annotated objects
% Created by Neeraj Kumar, please cite the following paper if you use this code-
% N. Kumar, R. Verma, S. Sharma, S. Bhargava, A. Vahadane and A. Sethi, 
% "A Dataset and a Technique for Generalized Nuclear Segmentation for
% Computational Pathology," in IEEE Transactions on Medical Imaging, 
% vol. 36, no. 7, pp. 1550-1560, July 2017
 
function [binary_mask,color_mask,JAI,DCoeff,JaccardIndex]=he_to_binary_mask(filename,folderXML,folderIMG,MS)
im_file=strcat(folderIMG,filename,'.tif');
 
xml_file=strcat(folderXML,filename,'.xml'); 
  
xDoc = xmlread(xml_file);
Regions=xDoc.getElementsByTagName('Region'); % get a list of all the region tags

MSLabel = bwlabel(MS);
MSLabelDice = MSLabel;
for regioni = 0:Regions.getLength-1
    Region=Regions.item(regioni);  % for each region tag
  
    %get a list of all the vertexes (which are in order)
    verticies=Region.getElementsByTagName('Vertex');
    xy{regioni+1}=zeros(verticies.getLength-1,2); %allocate space for them
    for vertexi = 0:verticies.getLength-1 %iterate through all verticies
        %get the x value of that vertex
        x=str2double(verticies.item(vertexi).getAttribute('X'));
        
        %get the y value of that vertex
        y=str2double(verticies.item(vertexi).getAttribute('Y'));
        xy{regioni+1}(vertexi+1,:)=[x,y]; % finally save them into the array
    end
    
end
im_info=imfinfo(im_file);
 
  
%nrow=im_info(s2).Height;
%ncol=im_info(s2).Width;
nrow=im_info.Height;
ncol=im_info.Width;

binary_mask=zeros(nrow,ncol); %pre-allocate a mask
color_mask = zeros(nrow,ncol,3);
%mask_final = [];
C=0;
U=0;
for zz=1:length(xy) %for each region
 %   fprintf('Processing object # %d \n',zz);
    smaller_x=xy{zz}(:,1); 
    smaller_y=xy{zz}(:,2);
    %make a mask and add it to the current mask
    %this addition makes it obvious when more than 1 layer overlap each
    %other, can be changed to simply an OR depending on application.
    polygon = poly2mask(smaller_x,smaller_y,nrow,ncol);
    
    
    binary_mask=binary_mask+zz*(1-min(1,binary_mask)).*polygon;%
    color_mask = color_mask + cat(3, rand*polygon, rand*polygon,rand*polygon);
    %binary mask for all objects
    %imshow(ditance_transform)
    %%Find segmented nuclei with major intersections
    nucleiIntersect = unique(MSLabel.*polygon);
    nucleiIntersect = nucleiIntersect(2:end);%2 because 0 its background
    JaccardIntersect =[];
    for ma = 1:length(nucleiIntersect)
        JaccardIntersect(ma)= jaccard(MSLabel==nucleiIntersect(ma)>0,polygon);
    end
    if ~isempty(JaccardIntersect)
    [Cupdate,idxMs] = max(JaccardIntersect);%Revisar Esto
    Uupdate = sumabs(or(MSLabel==nucleiIntersect(idxMs)>0,polygon));
    Cupdate = sumabs(and(MSLabel==nucleiIntersect(idxMs)>0,polygon));
    C = C+Cupdate;
    U = U+Uupdate;
    [idxF,idxR] = find(MSLabel==nucleiIntersect(idxMs));
    MSLabel(idxF,idxR)=0;
    end
    
    
    
end

    Uupdate = sumabs(MSLabel>0);
    U = U+Uupdate;
    
    JAI = C/U;
    
 DCoeff = dice(double(MS>0), double(binary_mask>0)); 
 JaccardIndex = jaccard(double(MS>0), double(binary_mask>0)); 
 
%figure;imshow(binary_mask)
%figure;imshow(color_mask)
