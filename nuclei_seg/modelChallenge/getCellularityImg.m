function feat = getCellularityImg( M )
    M = logical(M);        
    regionProperties = regionprops(M,'Centroid');
    nucleiCentroids = cat(1, regionProperties.Centroid);
    r=getSumNodeWeights(nucleiCentroids);
    feat=[length(nucleiCentroids) max(r) min(r) mean(r) median(r) std(r)];
end