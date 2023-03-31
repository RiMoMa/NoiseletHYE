


%tsne exploration
tiledlayout(2,2)
for K = [4,6,8,16]
fprintf('cluster %d\n',K)    
    
[idx,vocab] = kmeans(allFeatures,K,'distance',UsedDistance);
y = tsne(vocab,'distance','cityblock');
labelCluster = zeros(1,length(unique(idx)));
values = [];
for x=1:length(unique(idx))
    bool_id = idx==x;
    NucleiPercent = round(sum(allLabels(bool_id))/length(allLabels(bool_id))*100);
    
    fprintf('cluster %d, nuclei percent %d, nonnuclei %d \n',x, NucleiPercent,100-NucleiPercent);
    labelCluster(x) = NucleiPercent>50;
    if labelCluster(x) ==1
        values=[ values,NucleiPercent];
    else
        values=[ values,100-NucleiPercent];
    end
end

%figure;


%figure;

% Top plot
nexttile
%boxplot(values,labelCluster)
gscatter(y(:,1),y(:,2),labelCluster)
title([SelClass{1},' ',num2str(K),' ','Clusters'])

end
