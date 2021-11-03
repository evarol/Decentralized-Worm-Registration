function pc2fig(pc_cpd,int_matched,dotSize)

for t=1:length(pc_cpd)
    t
hold on
scatter3(pc_cpd{t}(:,1),pc_cpd{t}(:,2),pc_cpd{t}(:,3),dotSize,int_matched{t}(:,1),'filled');
axis equal;axis tight;grid on
drawnow
end
