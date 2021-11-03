function myplot3(X,args,size)
if nargin<3
    size=20;
end
plot3(X(:,1),X(:,2),X(:,3),args,'MarkerSize',size);
end