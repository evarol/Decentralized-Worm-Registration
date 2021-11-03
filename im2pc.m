function [pc,int]=im2pc(RFP,GFP,sigma)
I=RFP;
G=GFP;
foreground = and(imregionalmax(I),zscore(I,[],'all')>sigma);

L=bwlabeln(foreground);

r=regionprops3(L,I,'Centroid','MeanIntensity','WeightedCentroid','MaxIntensity');

pc=r.WeightedCentroid;
int(:,1)=r.MeanIntensity;

rg=regionprops3(L,G,'Centroid','MeanIntensity','WeightedCentroid','MaxIntensity');
int(:,2)=rg.MeanIntensity;

end




