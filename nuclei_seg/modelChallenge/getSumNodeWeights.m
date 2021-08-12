function vect=getSumNodeWeights( feature )
% GETSUMNODEWEIGHTS Computes the sum of the inverse of the distances 
% for each node in the graph
 
dist=pdist(feature);
dist(dist==0)=1;
dist=dist.^-1;
vect=sum(squareform(dist));   

end