function [cost,route]=error2path(ER,i0,j0)
cost=Inf;
t=0;
while cost==Inf
    t=t+1;
    G=zeros(size(ER));
    for i=1:size(ER,1)
        [~,idx]=sort(ER(i,:),'ascend','MissingPlacement','last');
        G(i,idx(1:t))=ER(i,idx(1:t));
        G(idx(1:t),i)=ER(idx(1:t),i);
        G(isnan(G))=0;
    end
    try;[cost,route]=dijkstra(G,i0,j0);end
end
end