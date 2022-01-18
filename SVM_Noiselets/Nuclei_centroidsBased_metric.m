function [Precision,Recall,Fscore]=Nuclei_centroidsBased_metric(groundT_mask, evaluated_mask) 

%Distance Thresholds for detected nuclei
tD = [2,4,8,16];
    
%determine centroids of ground truth
    groundT_centroids = regionprops(groundT_mask,'Centroid');% determine centroids of the ground truth mask
    N_nuclei = length(groundT_centroids); % number of annotated nuclei
    
    
    %%% make binary image a labeled image of the detected mask
    evaluated_mask_enumerated = bwlabeln (evaluated_mask); % label detected mask
    evaluated_centroids = regionprops(evaluated_mask_enumerated,'Centroid'); %find centroids of the detected mask 
    
    %% DISTANCES
    euclidean_distances = []; % array for computed distances
    
    FN= 0;
    coounter_distances = 1; %control counter for array distances 
    
    for n =1:N_nuclei %% detect and measure centroid distance for each annotatead nucleus
        
        Element_analyzed = (groundT_mask == n)*n ; % extract nucleus 
        Intersection_mask_b = and(Element_analyzed,evaluated_mask_enumerated); % determine intersected nuclei with the detected mask
        intersect_pixels = find(Intersection_mask_b==1); %find the pixel intersected idx  
        Intersect_elements = evaluated_mask_enumerated( intersect_pixels); %extract
        idx_element = unique(Intersect_elements); % determine nuclei that intersect the analyzed nucleus
        AreasIntersected = []; % array for area
        count=1; % counter for save the intersected areas
        
        for cA =idx_element'
            AreasIntersected(count) = sum(Intersect_elements(:) == cA); % determine the overlapping area for each nucleus
            count = count+1; %for next element
        end
        
        [~,Sel_idx]=max(AreasIntersected); % determine the biggest intersection
        Sel_idx = idx_element(Sel_idx); % determine the idx of the biggest intersect nucleus
        groundT_centroid_n = groundT_centroids(n).Centroid; %extract the centroid for distance compute, ground truth
        if ~isempty(Sel_idx)
            evaluated_centroid_sel = evaluated_centroids(Sel_idx).Centroid;
            evaluated_centroids(Sel_idx).Centroid = [NaN,NaN]; %delete element for future analysis
            nucleus_to_eliminate=evaluated_mask_enumerated==Sel_idx;%make a mask for eliminate nucleus for the original mask
            evaluated_mask_enumerated = evaluated_mask_enumerated-nucleus_to_eliminate; %final mask excluding for future analysis nucleus taked
            centroids_distance = pdist2(groundT_centroid_n,evaluated_centroid_sel,'euclidean');%compute distance between centroids
            euclidean_distances(coounter_distances) = centroids_distance; % array of distances
            coounter_distances = coounter_distances+1;
        else
            fprintf('none nucleus intersect ground truth\n')
            FN=FN+1;            
            
        end
        
    end
    try
    C = bsxfun(@lt,euclidean_distances',tD);%compute if is a TP for each treshold
    if all(size(C)==[1,4])
        C(2,:)=[0,0,0,0];
    end
    TP = sum(C); % sum of true positive nuclei
    catch
        TP=[NaN,NaN,NaN,NaN];
    end
    a = [evaluated_centroids.Centroid]; % for determine nucleus without intersection in the evaluated mask
    FP = sum(~ismissing(a))/2; %compute number of FP by counting NAN elements in a twoD array (->/2)
   
    %%%Final Performance metrics 
    Precision = TP./(TP+FP);
    Recall = TP./(FN+TP);
    Fscore = 2*(Precision.*Recall)./ (Precision+Recall);

end