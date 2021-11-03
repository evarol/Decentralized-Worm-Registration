function [atransform,distance,beta]=linhistmatch(a,b,nbins,type)
%takes as input two traces a (moving), b (reference) and outputs normalized
%trace atransform that has a similar histogram as b
%input : a - moving time series
%        b - reference time series
%        nbins - number of bins to discretize both time series into
%output: atransform - moved time series
%        distance - histogram distance between atransform and b

% discarding nans from time series
a_nan_idx=~isnan(a);
b_nan_idx=~isnan(b);
a=a(a_nan_idx);
b=b(b_nan_idx);
% Y=[fliplr(1-logspace(-6,0,nbins)) 1];
Y=linspace(0,1,nbins);
% discretizing time series using quantiles
abins=quantile(a,Y)'; 
bbins=quantile(b,Y)';
% weighted linear regression of the matching quantiles
if strcmpi(type,'non-negative')
    beta=lsqnonneg([abins ones(size(abins,1),1)],bbins); 
elseif strcmpi(type,'regular')
    beta=linsolve([abins ones(size(abins,1),1)],bbins);
% beta=pinv([abins ones(size(abins,1),1)])*bbins;
% x=[abins ones(size(abins,1),1)];
% beta=(x'*x + 0.01*eye(size(x'*x)))\(x'*bbins);
end
%transformed time series with nan's put back in

atransform=nan(size(a_nan_idx));
atransform(a_nan_idx)=a*beta(1) + beta(2);

%wasserstein distance computation
% atransformhat = histc(atransform,bbins)'; %discretizing the transformed time series
% distance=wdist(atransformhat,bhat,1); %wasserstein distance computation between atransform and b
distance=[]; %% turned this off since it seems to be not too useful
end
