condition = struct;
filename = "Dose_screen_fgfr1_pka.xls";

%% define conditions to test
% PKA up FGFR1 down
condition(1).nodes = [14,27];  
condition(1).inputvalues = [1,0];
condition(1).initial_attr=3 ; %Runx2+

% PKA up FGFR1 down
condition(2).nodes = [14,27]; 
condition(2).inputvalues = [1,0.2];
condition(2).initial_attr=3 ; %Runx2+

%% Alternative way to screen all possible values for PKA and FGFR1
% condi = 1;
% for pka= [0:0.1:1] %0:0.2:1
%     for fgfr=[0:0.1:1]
%         pka
%         fgfr
%         condition(condi).nodes = [14,27]; 
%         condition(condi).inputvalues = [pka,fgfr];
%         condition(condi).initial_attr=3 ; %Runx2+
%         condi = condi+1;
%     end
% end
%%


for c = 1 : length(condition)
    c
results_runs = Save_insilico_conditions_with_stats(condition(c).nodes,condition(c).inputvalues,condition(c).initial_attr,c,filename);
condition(c).result_transition = results_runs;
end


function [table2save]= Save_insilico_conditions_with_stats(nodes,inputvalues,initial_attr,sheetNumber,filename)
% run 3 times the same perturbed conditions in order to compute average and
% standard deviation in state transition.

runs = 3;
results_runs = zeros(runs,3);
for run = 1 : runs
    run
    [TransitionResult] = inSilicoConditions_return(2/3,'attractorAC_WT.mat',initial_attr,nodes,inputvalues);
    results_runs(run,:) = TransitionResult ;
end

results_runs(4,:) = mean(results_runs);
results_runs(5,:) = std(results_runs);

colNames = {'None','Sox9','Runx2'};
rowNames = {'run1','run2','run3','average','std'};
table2save = array2table(results_runs, 'VariableNames', colNames,'RowNames',rowNames);
writetable(table2save,filename,'WriteVariableNames',true,'WriteRowNames',true,'Sheet',sheetNumber,'Range','A2'); 

xlswrite(filename,[nodes;inputvalues],sheetNumber,'F2');




end