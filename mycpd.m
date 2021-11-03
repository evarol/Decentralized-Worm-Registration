function [moved,flow_field,cost]=mycpd(moving,fixed,sigma,iter,visual)
if nargin<5
    visual=0;
end
moving_orig=moving;
S=exp(-squareform(pdist(moving)).^2./2*sigma^2);
S=S./sum(S,2);
for t=1:iter
    D=pdist2(moving,fixed);
    [P,cost]=munkres(D);
    idx=find(any(P>0,2));
    vector = P(idx,:)*fixed-moving(idx,:);
    moving = moving + S(:,idx)*vector;
    if visual==1
        cla;hold on;
        myplot3(fixed,'b.',20);
        myplot3(moving,'r.',15);axis equal;axis tight;grid on;drawnow
    end
end
moved=moving;
flow_field = moved-moving_orig;


end