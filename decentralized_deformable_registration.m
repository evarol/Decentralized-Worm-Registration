function [pc_cpd,cpd_error,ER_cpd,cost_cpd,route_cpd,moved_cpd,flow_field_cpd]=decentralized_deformable_registration(pc_reg,demons_sigma,demons_iterations,cpd_subsampling_rate)
numT=length(pc_reg);
hold on
ER_cpd=nan(numT);
numsearch=numT^2;t=0;globalTic=tic;
bestsofar=Inf;
for i=1:numT
    for j=1:numT
        t=t+1;
        if or(and(rand(1)<=cpd_subsampling_rate*log(numT)/numT,i~=j),abs(i-j)<4)
            [moved_cpd{i,j},flow_field_cpd{i,j},ER_cpd(i,j)]=mycpd(pc_reg{i},pc_reg{j},demons_sigma,demons_iterations); %#ok<*AGROW>
            if ER_cpd(i,j)<bestsofar
                bestsofar=ER_cpd(i,j);
                cla;myplot3(moved_cpd{i,j},'b.');hold on;myplot3(pc_reg{j},'r.');axis equal;axis tight;grid on;drawnow
            end
        end
        clc
        fprintf(['Decentralized deformable registration - pairwise registration (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc(globalTic);
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
    end
end


[~,template_cpd]=min(nanmean(ER_cpd,2));
numsearch=numT;t=0;
for i=1:numT
    t=t+1;
    try
        if i==template_cpd
            cost_cpd{i,template_cpd}=0;
            route_cpd{i,template_cpd}=template_cpd;
        else
            [cost_cpd{i,template_cpd},route_cpd{i,template_cpd}]=error2path(ER_cpd,template_cpd,i);
        end
        pc_cpd{i}=pc_reg{i};
        for j=1:length(route_cpd{i,template_cpd})-1
            pc_cpd{i}=mycpd(pc_cpd{i},pc_reg{route_cpd{i,template_cpd}(j+1)},demons_sigma,demons_iterations);
        end
    catch
        warning('Route failed - skipping time point');
    end
    
    clc
    fprintf(['Decentralized deformable registration - Route centralization(' num2str(t) '/' num2str(numsearch) ')...\n']);
    fprintf(['\n' repmat('.',1,50) '\n\n'])
    for tt=1:round(t*50/(numsearch))
        fprintf('\b|\n');
    end
    TT=toc(globalTic);
    disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
end

numsearch=numT;t=0;globalTic=tic;
for i=1:numT
    t=t+1;
    [~,cpd_error(i,1)]=munkres(pdist2(pc_cpd{template_cpd},pc_cpd{i}));
    clc
    fprintf(['Decentralized deformable registration - Error computing (' num2str(t) '/' num2str(numsearch) ')...\n']);
    fprintf(['\n' repmat('.',1,50) '\n\n'])
    for tt=1:round(t*50/(numsearch))
        fprintf('\b|\n');
    end
    TT=toc(globalTic);
    disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
end