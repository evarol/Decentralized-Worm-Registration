function cs=pc2cs(pc,margin)
%find minimal set such that points are at most margin away from any coreset
%element

[~,cliques]=maximal_clique(squareform(pdist(pc))<margin,size(pc,1),1);

cs=[];
for i=1:length(cliques)
    if ~isempty(cliques{i})
      cs=[cs;mean(pc(cliques{i},:),1)];
    end
end

end