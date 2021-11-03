function out=ts2gp(ts,L,sigma)

moving=(1:length(ts))';
fixed=(1:length(ts))';

K1=exp(-squareform(pdist(moving).^2/2*L^2)) + eye(size(moving,1))*sigma.^2;
K2=exp(-pdist2(moving,fixed).^2/2*L^2);

out = (K2'/K1)*ts;


end

