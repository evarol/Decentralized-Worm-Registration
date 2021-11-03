function out=pc2gp(moving,fixed,moving_int,L,sigma)

K1=exp(-squareform(pdist(moving,'mahalanobis').^2/2*L^2)) + eye(size(moving,1))*sigma.^2;
K2=exp(-pdist2(moving,fixed,'mahalanobis').^2/2*L^2);

out = (K2'/K1)*moving_int;

end