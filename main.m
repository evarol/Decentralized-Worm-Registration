clearvars -except b pc cs cs2 R T ER pc_reg unreg_error reg_error tform ER_cpd cpd_error int pc_cpd int_gp moved_cpd ER_cpd flow_field_cpd route route_cpd cost_cpd cost int_matched
close all
clc
minmax = @(x)((x-nanmin(x(:)))./nanmax(x(:)-nanmin(x(:))));
addpath(genpath('/Users/erdem/Dropbox/matlab_toolboxes/export-fig'));
startTic=tic;


%% pixel resampling parameters
% pixel_dimensions=[0.27 0.27 1.5];
pixel_dimensions=[0.4 0.4 1.5];
upsample_factor=1;

%% Detection parameters
sigma=16;


%% Rigid registration parameters
rigid_subsampling_rate=1;
cpd_subsampling_rate=1;

%% Deformable registration parameters
demons_sigma=0.05;
demons_iterations=10;

%% Signal standardization parameters
% gp_L=0.1;
% gp_sigma=0.1;

gp_L=1;
gp_sigma=1;
%% Visualization parameters
rasterSize=0.1;
dotSize=10;


if ~exist('b');
    b=h5read('./data/data_1.h5','/data');
    tmp=b(:,:,:,1,:);
    b(:,:,:,1,:)=b(:,:,:,2,:);
    b(:,:,:,2,:)=tmp;clear tmp
    b=b(:,:,:,:,1:100);
end
numT=size(b,5);


if ~exist('pc')
    numsearch=size(b,5);globalTic=tic;
    for t=1:numT
        [pc{t},int{t}]=im2pc(imresize3(double(b(:,:,:,1,t)),size(b(:,:,:,1,1)).*pixel_dimensions*upsample_factor),imresize3(double(b(:,:,:,2,t)),size(b(:,:,:,2,1)).*pixel_dimensions*upsample_factor),sigma); %% maybe first resample to um before getting detections
        pc{t}=pc{t}./upsample_factor;
        clc
        fprintf(['Extracting point clouds (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc(globalTic);
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
        
    end
end

if ~exist('R');
    hold on
    ER=nan(numT);   
    numsearch=numT*(numT-1)/2;t=0;globalTic=tic;
    bestsofar=Inf;
    for i=1:numT-1
        for j=i+1:numT
            t=t+1;
            if or(rand(1)<=rigid_subsampling_rate*log(numT)/numT,j<i+4)
                [rot, trans, err, ~] = icp2(pc{i}',pc{j}');
                R{i,j}=rot';
                T{i,j}=trans';
                R{j,i}=R{i,j}';
                T{j,i}=-T{i,j}*R{i,j}';
                ER(i,j)=err(end);ER(j,i)=err(end);
                if ER(i,j)<bestsofar
                    bestsofar=ER(i,j);
                    cla;myplot3(pc{j}*R{i,j}+T{i,j},'b.');hold on;myplot3(pc{i},'r.');axis equal;axis tight;grid on;drawnow
                end
            end
            clc
            fprintf(['Decentralized registration (' num2str(t) '/' num2str(numsearch) ')...\n']);
            fprintf(['\n' repmat('.',1,50) '\n\n'])
            for tt=1:round(t*50/(numsearch))
                fprintf('\b|\n');
            end
            TT=toc(globalTic);
            disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
        end
    end
    
end


[~,template]=min(nanmean(ER,2));
if ~exist('pc_reg')
    for i=1:numT
        i
        try
            if i==template
                cost{i,template}=0;
                route{i,template}=template;
            else
                [cost{i,template},route{i,template}]=error2path(ER,template,i);
            end
            pc_reg{i}=pc{i};
            for t=1:length(route{i,template})-1
                pc_reg{i}=pc_reg{i}*R{route{i,template}(t+1),route{i,template}(t)} + T{route{i,template}(t+1),route{i,template}(t)};
            end
        end
    end
    
    
    for i=1:numT
        i
        [~,unreg_error(i,1)]=munkres(pdist2(pc{template},pc{i}));
        [~,reg_error(i,1)]=munkres(pdist2(pc_reg{template},pc_reg{i}));
    end
end

if ~exist('moved_cpd');
    hold on
    ER_cpd=nan(numT);   
    numsearch=numT^2;t=0;globalTic=tic;
    bestsofar=Inf;
    for i=1:numT
        for j=1:numT
            t=t+1;
            if or(and(rand(1)<=cpd_subsampling_rate*log(numT)/numT,i~=j),abs(i-j)<4)
                [moved_cpd{i,j},flow_field_cpd{i,j},ER_cpd(i,j)]=mycpd(pc_reg{i},pc_reg{j},demons_sigma,demons_iterations);
                if ER_cpd(i,j)<bestsofar
                    bestsofar=ER_cpd(i,j);
                    cla;myplot3(moved_cpd{i,j},'b.');hold on;myplot3(pc_reg{j},'r.');axis equal;axis tight;grid on;drawnow
                end
            end
            clc
            fprintf(['Decentralized deformable registration (' num2str(t) '/' num2str(numsearch) ')...\n']);
            fprintf(['\n' repmat('.',1,50) '\n\n'])
            for tt=1:round(t*50/(numsearch))
                fprintf('\b|\n');
            end
            TT=toc(globalTic);
            disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
        end
    end
    
end

[~,template_cpd]=min(nanmean(ER_cpd,2));
if ~exist('pc_cpd')
    for i=1:numT
        i
        try
            if i==template_cpd
                cost_cpd{i,template_cpd}=0;
                route_cpd{i,template_cpd}=template_cpd;
            else
                [cost_cpd{i,template_cpd},route_cpd{i,template_cpd}]=error2path(ER_cpd,template_cpd,i);
            end
            pc_cpd{i}=pc_reg{i};
            for t=1:length(route_cpd{i,template_cpd})-1
                pc_cpd{i}=mycpd(pc_cpd{i},pc_reg{route_cpd{i,template_cpd}(t+1)},demons_sigma,demons_iterations);
            end
        end
    end
    
    
    for i=1:numT
        i
        [~,cpd_error(i,1)]=munkres(pdist2(pc_cpd{template_cpd},pc_cpd{i}));
    end
end
% 
% if ~exist('pc_cpd');
%     numsearch=numT;t=0;globalTic=tic;
%     for i=1:numT
%         t=t+1;
%         moving=pc_reg{i};
%         fixed=pc_reg{template};
%         pc_cpd{i}=mycpd(moving,fixed,demons_sigma,demons_iterations);
%          clc
%             fprintf(['Centralized deformable registration (' num2str(t) '/' num2str(numsearch) ')...\n']);
%             fprintf(['\n' repmat('.',1,50) '\n\n'])
%             for tt=1:round(t*50/(numsearch))
%                 fprintf('\b|\n');
%             end
%             TT=toc(globalTic);
%             disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
%     end
%     
%     for i=1:numT
%         i
%         [~,cpd_error(i,1)]=munkres(pdist2(pc_cpd{template},pc_cpd{i}));
%     end
% end


if ~exist('int_gp');
    for i=1:numT
        int_matched{i}=linhistmatch(int{i},int{template},20,'regular');
        moving=[pc_cpd{i} int_matched{i}];
        fixed=[pc_cpd{template} int{template}];
        moving_int=int_matched{i};
        int_gp{i}=pc2gp(moving,fixed,moving_int,gp_L,gp_sigma);
        i
    end
end

for t=1:length(int_gp)
X(:,t,1)=int_gp{t}(:,1);
X(:,t,2)=int_gp{t}(:,2);
end
stopToc=toc(startTic);
disp(['Time taken: ' num2str(stopToc/60) ' minutes.']);

figure('units','normalized','outerposition',[0 0 1 1])
for t=1:length(pc_cpd)
    t
hold on
scatter3(pc_cpd{t}(:,1),pc_cpd{t}(:,2),pc_cpd{t}(:,3),dotSize,int_matched{t}(:,1),'filled');
axis equal;axis tight;grid on
drawnow
end

rigid_route=zeros(length(pc));
cpd_route=zeros(length(pc));
for t=1:size(route,1)
for i=1:length(route{t,template})-1
rigid_route(route{t,template}(i),route{t,template}(i+1))=1;
end
end


for t=1:size(route_cpd,1)
for i=1:length(route_cpd{t,template_cpd})-1
cpd_route(route_cpd{t,template_cpd}(i),route_cpd{t,template_cpd}(i+1))=1;
end
end

g_rigid=digraph(rigid_route);
g_cpd=digraph(cpd_route);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
h=plot(g_rigid,'Layout','force','MarkerSize',10,'NodeColor','b','EdgeColor','k','NodeFontSize',14);
hold on;plot(h.XData(template),h.YData(template),'g.','MarkerSize',40)
legend('Rigid routes','template')
subplot(1,2,2)
% plot(g_cpd,'XData',h.XData,'YData',h.YData,'MarkerSize',10,'NodeColor','r','EdgeColor','k','NodeFontSize',14);
h=plot(g_cpd,'Layout','force','MarkerSize',10,'NodeColor','r','EdgeColor','k','NodeFontSize',14);
hold on;plot(h.XData(template_cpd),h.YData(template_cpd),'g.','MarkerSize',40)
legend('Deformable routes','template')


figure('units','normalized','outerposition',[0 0 1 1])
for i=1:numT
    subplot(2,3,1);myplot3(pc{i},'b.');set(gca,'xlim',[0 500*pixel_dimensions(1)],'ylim',[0 500*pixel_dimensions(2)],'zlim',[0 100*pixel_dimensions(3)]);grid on;title(['Unregistered - Frame ' num2str(i)]);set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');drawnow
    subplot(2,3,2);myplot3(pc_reg{i},'r.');set(gca,'xlim',[0 500*pixel_dimensions(1)],'ylim',[0 500*pixel_dimensions(2)],'zlim',[0 100*pixel_dimensions(3)]);grid on;title(['Coarse Rigid Decentralized Registered to frame ' num2str(template) ' - Frame ' num2str(i)]);set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');drawnow
    subplot(2,3,3);myplot3(pc_cpd{i},'r.');set(gca,'xlim',[0 500*pixel_dimensions(1)],'ylim',[0 500*pixel_dimensions(2)],'zlim',[0 100*pixel_dimensions(3)]);grid on;title(['Deformable Registered to frame ' num2str(template) ' - Frame ' num2str(i)]);set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');drawnow
    subplot(2,3,[4:6]);plot(1:i,unreg_error(1:i),'b',1:i,reg_error(1:i),'r',1:i,cpd_error(1:i),'g','LineWidth',2);legend('Unregistered error','Rigid registered error','Deformable registered error');set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');grid on;drawnow
    set(gcf,'Color','w');
%     export_fig(['./figs/frame_' num2str(i) '.png']);
end

figure;
subplot(4,1,1)
hold on
for t=1:length(pc_reg)
    scatter(pc{t}(:,1)+t*512*pixel_dimensions(2),pc{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,2),'filled');
    drawnow
end
colorbar
title('Unregistered Postural Raster Plot - GFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,2)
hold on
for t=1:length(pc_reg)
    scatter(pc_reg{t}(:,1)+t*512*pixel_dimensions(2),pc_reg{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,2),'filled');
    drawnow
end
colorbar
title('Rigid registered Postural Raster Plot - GFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,3)
hold on
for t=1:length(pc_reg)
    scatter(pc_cpd{t}(:,1)+t*512*pixel_dimensions(2),pc_cpd{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,2),'filled');
    drawnow
end
colorbar
title('Deformable registered Postural Raster Plot - GFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,4)
hold on
for t=1:length(pc_reg)
    scatter(pc_cpd{template}(:,1)+t*512*pixel_dimensions(2),pc_cpd{template}(:,2),size(pc{template},1)*rasterSize,int_gp{t}(:,2),'filled');
    drawnow
end
colorbar
title('GP interpolated GFP intensities on template ');
xlabel('Horizontal space (um) x time(s)');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');


figure;
subplot(4,1,1)
hold on
for t=1:length(pc_reg)
    scatter(pc{t}(:,1)+t*512*pixel_dimensions(2),pc{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,1),'filled');
    drawnow
end
colorbar
title('Unregistered Postural Raster Plot - RFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,2)
hold on
for t=1:length(pc_reg)
    scatter(pc_reg{t}(:,1)+t*512*pixel_dimensions(2),pc_reg{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,1),'filled');
    drawnow
end
colorbar
title('Rigid registered Postural Raster Plot - RFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,3)
hold on
for t=1:length(pc_reg)
    scatter(pc_cpd{t}(:,1)+t*512*pixel_dimensions(2),pc_cpd{t}(:,2),size(pc{t},1)*rasterSize,int{t}(:,1),'filled');
    drawnow
end
colorbar
title('Deformable registered Postural Raster Plot - RFP intensities');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
subplot(4,1,4)
hold on
for t=1:length(pc_reg)
    scatter(pc_cpd{template}(:,1)+t*512*pixel_dimensions(2),pc_cpd{template}(:,2),size(pc{template},1)*rasterSize,int_gp{t}(:,1),'filled');
    drawnow
end
colorbar
title('GP interpolated RFP intensities on template ');
xlabel('Horizontal space (um) x time(s)');
ylabel('Vertical space (um)');
set(gca,'FontWeight','bold','FontSize',20,'TickLength',[0 0]);set(gcf,'Color','w');
