function [pc_reg,reg_error,unreg_error,ER,cost,route,template]=decentralized_rigid_registration(pc,rigid_subsampling_rate)
numT=length(pc);
hold on
ER=nan(numT);
numsearch=numT*(numT-1)/2;t=0;globalTic=tic;
bestsofar=Inf;
for i=1:numT-1
    for j=i+1:numT
        t=t+1;
        if or(rand(1)<=rigid_subsampling_rate*log(numT)/numT,j<i+4)
            [rot, trans, err, ~] = icp2(pc{i}',pc{j}');
            R{i,j}=rot'; %#ok<*AGROW>
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
        fprintf(['Decentralized registration - pairwise registration (' num2str(t) '/' num2str(numsearch) ')...\n']);
        fprintf(['\n' repmat('.',1,50) '\n\n'])
        for tt=1:round(t*50/(numsearch))
            fprintf('\b|\n');
        end
        TT=toc(globalTic);
        disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
    end
end



[~,template]=min(nanmean(ER,2));
numsearch=numT;t=0;globalTic=tic;
for i=1:numT
    t=t+1;
    try
        if i==template
            cost{i,template}=0;
            route{i,template}=template;
        else
            [cost{i,template},route{i,template}]=error2path(ER,template,i);
        end
        pc_reg{i}=pc{i};
        for j=1:length(route{i,template})-1
            pc_reg{i}=pc_reg{i}*R{route{i,template}(j+1),route{i,template}(j)} + T{route{i,template}(j+1),route{i,template}(j)};
        end
    catch
        warning('Route failed - skipping time point');
    end
    clc
    fprintf(['Decentralized registration - Route centralization (' num2str(t) '/' num2str(numsearch) ')...\n']);
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
    [~,unreg_error(i,1)]=munkres(pdist2(pc{template},pc{i}));
    [~,reg_error(i,1)]=munkres(pdist2(pc_reg{template},pc_reg{i}));
    clc
    fprintf(['Decentralized registration - Error computing (' num2str(t) '/' num2str(numsearch) ')...\n']);
    fprintf(['\n' repmat('.',1,50) '\n\n'])
    for tt=1:round(t*50/(numsearch))
        fprintf('\b|\n');
    end
    TT=toc(globalTic);
    disp(['Time elapsed (minutes): ' num2str(TT/60) ' Time remaining (minutes): ' num2str((numsearch-t)*(TT/t)*(1/60)) ' Est. Total (minutes): ' num2str(TT/60 + (numsearch-t)*(TT/t)*(1/60))]);
end
end
