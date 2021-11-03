function [lb,bestset]=pc2pc(pc1,pc2,margin,radius,maxdepth)
addpath(genpath('/Users/erdem/Documents/Github/Maximal-rotation-set'));
D1=squareform(pdist(pc1));
D2=squareform(pdist(pc2));

D1(D1>radius)=Inf;
D2(D2>radius)=Inf;

lb=maximal_rotation_set_DFS(D1,D2,margin,maxdepth);

[i0,j0]=find(lb==max(lb(:)));i0=i0(1);j0=j0(1);
bestset=maximal_rotation_set_DFS_decode(D1,D2,margin,maxdepth,i0,j0);

end