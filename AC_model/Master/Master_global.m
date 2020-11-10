%% Master script combining AtoZ analysis from publication %%

% rng('shuffle')
global n
n = 60;

global tolerance
tolerance = 1e-3;

%% Compute WT attractors
    disp('Compute WT attractors...')
    nstates = 10;
    saturation = 2/3;
    input = [];
    inputindex = [];
    [d, sizes, da, percent, attractors,attractors_read,T_sizes,T_percent,Attr_profiles] = Attractor_AC(nstates, saturation,input, inputindex);

    string = 'attractorAC_WT'; 
    save([string,'_readerFriendly.mat'], 'T_sizes','T_percent','Attr_profiles');

%% Compute constrained attractors (growth factor profiles)
    disp('Compute constrained attractors...')
    
    % Healthy like environment:
    inputindex = [1,3,5,11,23,33,51]; %growth factor indexes
    input = [0,1,0.2,1,0.7,0,0]; % imposed profile
    [~,sizes_healthy, ~, percent_healthy, attractors_healthy,attractors_read_healthy,T_sizes_healthy,T_percent_healthy,Attr_profiles_healthy] = Attractor_AC(nstates, saturation,input, inputindex);

    string = 'attractorAC_HealthyCues'; 
    save([string,'.mat'], 'T_sizes_healthy','T_percent_healthy','Attr_profiles_healthy');

    % Hypertrophy like environment:
    inputindex = [1,3,5,11,23,33,51]; %growth factor indexes
    input = [1,0,1,0,0,1,1]; % imposed profile
    [~,sizes_oa, ~, percent_oa, attractors_oa,attractors_read_oa,T_sizes_oa,T_percent_oa,Attr_profiles_oa] = Attractor_AC(nstates, saturation,input, inputindex);

    string = 'attractorAC_OAcues'; 
    save([string,'.mat'], 'T_sizes_oa','T_percent_oa','Attr_profiles_oa');
   
%% Single perturbation
    disp('Screen single node perturbations...')
    [perresult, stabattr, nb_transition,storeNewattr,transition,twicepresult,distribution] = StateStab_perturbation(saturation,attractors);
    % save variables in file.mat:
    string = ['SingleNodePertub' '_sat' num2str(saturation) '.mat'];
    save(string, 'stabattr', 'perresult', 'trlik', 'trlik2','saturation','storeNewattr','transition','distribution')

%% Test any in silico condition and compute descriptive statistics
% NB: not entirely robust for different WT attractors, adaptation required
    disp('Specific case study - in silico condition...')
    % 1/ define conditions to test: e.g all possible values of pka&fgfr1 combi
    condition = struct;
    condi = 1;
    for pka= 0:0.1:1 %0:0.2:1
        for fgfr= 0:0.1:1
            pka;
            fgfr;
            condition(condi).nodes = [14,27]; 
            condition(condi).inputvalues = [pka,fgfr];
            condition(condi).initial_attr=2 ; %Runx2+
            condi = condi+1;
        end
    end
    
    % 2/ compute transitions for each condition of interest
    filename = "In_silico_conditions.xls";
    for c = 1 : length(condition)
    results_runs = Save_insilico_conditions_with_stats(attractors,condition(c).nodes,condition(c).inputvalues,condition(c).initial_attr,c,filename,saturation);
    condition(c).result_transition = results_runs;
    end

%% Screen all pairwise combinations
%possibility to run this part in parallel (for HPC)
    disp('Screen pairwise perturbations...')
    disp('length of attractrors is:')
    disp(length(attractors));
    
    currentFolder = pwd;
    path_dir =  'Screening_results';
    mkdir Screening_results;
    
    C=combnk(1:n,2); % all combinations of 2 among 60. 
    %disp(['matlab number of cores : '   num2str(feature('numCores'))])

    for comb=1:size(C,1)
        Ccomb=C(comb,:);
        nod1=Ccomb(1);
        nod2= Ccomb(2);
        Bi_perturbations_noNone(saturation,1,attractors,nod1,nod2,path_dir);   
    end
    
% Summarize pairwise screening results and analyse:
    disp('Summarize and analyze pairwise combination screening results...')
    [Summary_screening,MyCombi] = screening_summary(C,path_dir);
    h = plot_PCA_pertub_screening(Summary_screening,MyCombi)
    
%% Select perturbations based on transition results    
% Used to select conditions leading the most to Sox9+ when starting from Runx2+
%Warning : This part is neither robust nor optimized so it requires manual
%adapatations to be apllied to other models or cases. 
    disp('Select perturbations of interest based on manually defined transition thresholds...')
    
    for t=70:10:100    %define threshold range
    %% parameters
    percentSox9= t; %threshold for % of transitions to Sox9+ 
    percentRunx2 = 100-percentSox9; %threshold for % of transitions to Runx2+ 
    %%
    init_attr=["Sox9" "Runx2"]; %["Sox9","Runx2"];
    SelectedCombi=select_perturbations(C,init_attr,percentSox9,percentRunx2,path_dir);
    save(['Combi2test' num2str(percentSox9) '_' num2str(percentRunx2) '.mat'],'SelectedCombi')
    end
    
    disp('write conditions selected based on defined threshold...')
    Report_selected_conditions() %
    disp('Output in excel file')
    disp('Done!')
    
%% Tool functions

% Monte Carlo utils functions:
function [d, sizes, da, percent, attractors,attractors_read,T_sizes,T_percent,Attr_profiles] = Attractor_AC(nstates, saturation,input, inputindex)
% State space screening with Monte Carlo analysis to find attractors
% 1 argument + 3 optional arguments
% The script requires: 
%   - File with equations for asynchronous updating: zupdateas_AC.m
%   - List of biological entity names and indexes: components.mat

c = 1; % activity = c*PTM*transcription/translation/degradation
% set optional parameters
    if nargin == 1
        input = [];
        inputindex = [];
        saturation = 2/3;
    end
    if nargin == 2
            input = [];
        inputindex = [];
    end

    if ~isnumeric(nstates)
       nstates = str2double(nstates);
    end

global n
global tolerance
% saturation = 1;
global s
s = zeros(1,5);
s(1) = 1;
for i = 2:5
s(i) = 2* saturation/i;
end

% define idex of nodes regulated at fast or slow level only: 
    slowidx = [3,5,11,13,16,20,21,23,24,25,44,46,47,48,49,51,59]; 
    fastidx = [2, 7, 14, 18, 19, 22, 26, 29, 30, 35, 36, 40, 41, 45, 50, 54, 55, 56];

d = zeros(1,nstates); % distance vectore: number of loops to reach steady state

k = 1;  % k is the attractor index (increment)
ptr = 1; % counts the total amount of initial states currently assessed

% initialisize some variables
    attractors = {ones(35,1)}; % will be used to store the attractors
    attrcnt = 0;
    aSox9 = 0;
    aRunx2 = 0;
    aMMP13= 0;
    
wb = waitbar(0,'computing attractors and distances...');

while ptr <= nstates 
    ptr;
    waitbar(ptr/nstates,wb);
    
    % generate random new initial state:
    newstate = zeros(3,n); 
    newstate(1,:) = rand(n,1); %first row= global activity
    newstate(1,inputindex) = input; % keep inputnodes fixed
    newstate(2,:) = sqrt(newstate(1,:)); %second row= fast/protein level
    newstate(3,:) = newstate(2,:); %third row= slow/gene level
    newstate(2,slowidx) = 1; % full control of fast controlled nodes in fast variable, same for slow
    newstate(3,slowidx) = newstate(1,slowidx);
    newstate(3,fastidx) = 1;
    newstate(2,fastidx) = newstate(1,fastidx);

    attr = 0;
    cnt = 1; % indicates the number of the state in the new sequence (reinitialises after every new random init state)
 while (attr == 0)
    attr = 1;
    for l = randperm(n)
        attrf = 0;
        while (attrf == 0)
            attrf = 1;
            for j = randperm(n)
                changednode = zupdateas_AC(newstate,j,'fast'); 
                
            
            if (abs(newstate(2,j)-changednode) > tolerance) && ~ismember(j,inputindex) % Check if new state is different from previous one
               newstate(2,j) = changednode;
               newstate(1,j) =  min(1,c*newstate(2,j)*newstate(3,j)); % min(newstate(2,j), newstate(3,j));
               attrf = 0;
               newstate(1,inputindex) = input; % keep inputnodes fixed
               break;
            end
            end
        end
            changednode = zupdateas_AC(newstate,l,'slow');
            
         if (abs(newstate(3,l)-changednode) > tolerance) && ~ismember(l,inputindex) % Check if new state is different from previous one
           newstate(3,l) = changednode;
           newstate(1,l) = min(1,c*newstate(2,l)*newstate(3,l));
           attr = 0;
           cnt = cnt + 1;
           newstate(1,inputindex) = input; % keep inputnodes fixed          
           break;
         end
    end     
 end

   % 2 options: we encountered attractor before, in which case we need to
   % know which k it corresponds to or we found a new attractor.
   
   if checkstates(newstate,attractors,k-1) % We encountered attractor before
      [~,old] = checkstates(newstate,attractors,k); 
            d(ptr) = cnt;   % record the distance to the attractor
            attrcnt(old) = attrcnt(old) + 1;

   else % We found a new attractor
       d(ptr) = cnt; 
       attractors(k) = {newstate};
       attrcnt(k) = 1;
       % store value of some common markers for the new attractor : 
       aRunx2(k)=newstate(1,9);
       aSox9(k)=newstate(1,10);
       aMMP13(k)=newstate(1,24);
       k = k + 1;
   end
   ptr = ptr + 1;
end

close(wb)

%% Determine size basins and values of Sox9/Runx2

k = k-1; % Since we never found an attractor to match the last k

sizes = zeros(6,k);
sizes(1,1:k) = 1:k; % numbering the attractors
sizes(3,:) = attrcnt;
sizes(4,:) = aSox9; % Store the amount of init states reaching attractors that are Sox9 positive
sizes(5,:) = aRunx2; 
sizes(6,:) = aMMP13;

% percentages in Runx2, Sox9, both, none
percent = cell(4,2);
percent(1,1) = {'Runx2'};
percent(2,1) = {'Sox9'};
percent(3,1) = {'both'};
percent(4,1) = {'none'};
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;

% attractors are positive for a marker when its value is > a threshold.    
presence_threshold = 0.05;
for i = 1:k
    if sizes(4,i) > presence_threshold
        if sizes(5,i) > presence_threshold
            temp3 = temp3 + sizes(3,i)/sum(sizes(3,1:k));
        else
            temp2 = temp2 + sizes(3,i)/sum(sizes(3,1:k));
        end
    else
        if sizes(5,i) > presence_threshold
            temp1 = temp1 + sizes(3,i)/sum(sizes(3,1:k));
        else
            temp4 = temp4 + sizes(3,i)/sum(sizes(3,1:k));
        end
    end
end
    percent(1,2) = {temp1};
    percent(2,2) = {temp2};
    percent(3,2) = {temp3};
    percent(4,2) = {temp4}; 
    
% Average distance to attractor: algorithm does not count states that start on
% an attractor, but is only insignificant fraction of state space.
da = zeros(size(d));
    da(1) = d(1);
for i= 1:length(d) - 1
    if d(i) <= d(i+1) 
        da(i) = d(i+1);
    end
end
da(da < 1) = [];
da = mean(da);


%% saving outputs in a readable format
%string = 'attractorAC_WT'; 

[attractors_read,T_sizes,T_percent,Attr_profiles] = reader_friendly_output(sizes, percent, attractors);

%save([string,'.mat'], 'd', 'sizes', 'da', 'percent', 'attractors','attractors_read');
%save([string,'_readerFriendly.mat'], 'T_sizes','T_percent','Attr_profiles');

end

function [boolean, number] = checkstates(newstate,state,maxidx)
% function to check whether a state has already been described
global tolerance

boolean = 0;
number = 0;

if maxidx ~= 0
    for i = 1:maxidx
        if  abs(newstate-state{i}) < 100*tolerance % isequal(newstate, state{i});
            boolean = 1;
            number = i;
            
            break;
        end
    end
end
end

function [attractors_read,T_sizes,T_percent,Profiles] = reader_friendly_output(sizes, percent, attractors)
    
    attractors_read=attractors;
    % formatting for easy reading:
        load('components.mat','componentnames')
        for a = 1:length(attractors_read)
            attractors_read{a}=[componentnames{1,:}; attractors_read{a}];
        end
% Table of attractor sizes:
    Attr_index = sizes(1,:)';
    Canalization = sizes(3,:)';
    ValueSox9 = sizes(4,:)';
    ValueRunx2 = sizes(5,:)';
    ValueMMP13 = sizes(6,:)';
    T_sizes = table(Attr_index, Canalization,ValueSox9,ValueRunx2,ValueMMP13);
    
% Table of percentage
    AttractorType = percent(:,1);
    PercentageInitialisations = percent(:,2);
    T_percent = table(AttractorType,PercentageInitialisations);
    
% Attractor profiles

Profiles = struct();
    for a = 1:length(attractors) 
        NodeNames = attractors_read{a}(1,:)';
        Global_activity = attractors{a}(1,:)';    
        Protein_activation = attractors{a}(2,:)';
        Gene_expression = attractors{a}(3,:)';
        Profiles(a).Attr_index= table(NodeNames,Global_activity,Protein_activation,Gene_expression);
    end

end

% Single node perubrations associated functions:
function [perresult, stabattr, nb_transition,storeNewattr,transition,twicepresult,distribution] = StateStab_perturbation(saturation,attractors)
% StateStab_perturbation: Calculates the stability of a states by perturbing network nodes one by
% one

%set parameters:
runs = 100; %number of repetition of the perturbation from init-state
fixedtime = 1000; %duration of perturbation (in time steps)

if ~isnumeric(saturation)
   saturation = str2num(saturation);
end

global n
   
global s
s = zeros(1,5);
s(1) = 1;
for i = 2:5
s(i) = 2* saturation/i;
end

perresult = cell(length(attractors),n,2); 
storeNewattr= cell(length(attractors),n,2);
twicepresult = cell(length(attractors),n,2); 

tic

for i = 1:length(attractors)
    disp(['attractor ' num2str(i) ' started...'])
    attr = attractors{i}; 
    
    for j = [1:n] 
        
        disp(['node ' num2str(j) ' ...'])
        if (attr(1,j) ~= 1) && (attr(1,j) ~= 0) % 2 possible flip states: 0 and 1: first 1 % rest algoritme uitwerken
        for k = 1:2
            per = attr;
            per(1,j) = k-1;% create perturbed state
            if ~isequal(attr,zupdates_AC(per))
                [attrcnt,~,stateruns] = ssattr(per, runs, attractors, fixedtime,j);
                perresult(i,j,k) = {attrcnt/2}; % store results of perturbation
                twicepresult(i,j,k) = {attrcnt};
                storeNewattr(i,j,k)={stateruns};  
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                twicepresult(i,j,k) = {attrcnt};
                attrcnt(i) = runs/2; % we want to give this the same weight: since every node has weight one, and we perform 2 simulations here, we divide by 2
                perresult(i,i,k) = {attrcnt};
                storeNewattr(i,j,k)={'no transition'};  
            end
            
        end
        else % we explore 2 values  only if the initial one is intermediate
        per = attr;
        per(1,j) = ~per(1,j);% create perturbed state
            if ~isequal(attr,zupdates_AC(per))
                [attrcnt,~,stateruns] = ssattr(per, runs, attractors, fixedtime,j);
                twicepresult(i,j,int8(-attr(1,j)+2)) = {attrcnt};
                perresult(i,j,int8(-attr(1,j)+2)) = {attrcnt}; % store results of perturbation
                storeNewattr(i,j,int8(-attr(1,j)+2))={stateruns};
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                attrcnt(i) = runs; 
                twicepresult(i,j,int8(-attr(1,j)+2)) = {attrcnt};
                perresult(i,i,int8(-attr(1,j)+2)) = {attrcnt};
                storeNewattr(i,j,int8(-attr(1,j)+2))={stateruns};
            end
        
        end
    
    end
    
end
% 

nattr = length(attractors); % determine the number of attractors
stabattr = cell(3,nattr,2); % stores how many flips have the potential to change attractors, second row contains overall prob, third row probability per state
nb_transition = zeros(nattr,nattr,2); % contains likelihood of transition (going from attr to another attr)
for c = 1:2
for a = 1:nattr
     flip = 0; % contains the number of flipped genes
     perflip = zeros(1, n); % contains the likelihood of flipping per gene
     for b = 1:n
        if ~isempty(twicepresult{a,b,c})
        x = twicepresult{a,b,c};
        if (x(a) ~= runs)
            flip = flip + 1;
            perflip(b) = (sum(x) - x(a))/runs; % percentage of states flipping
            nb_transition(a,1:nattr,c) = nb_transition(a,1:nattr,c) + x(1:nattr); % add number of transitions to transition matrix
        else
            nb_transition(a,a,c)= nb_transition(a,a,c) + runs; 
        end
        end
     end

     stabattr(1,a,c) = {flip};
     stabattr(3,a,c) = {perflip};
     stabattr(2,a,c) = {sum(perflip)/(flip*runs)}; % average chance of flipping per node that can potentially switch states
end
end
toc
tr_likood = mean(nb_transition,3)*2;

%calculation of matrix of transition:
transition=tr_likood./sum(tr_likood,2); 
distribution=stationary(transition);

end

function [attrcnt,newattr,stateruns] = ssattr(initstate,runs,attractors, fixedtime, fixednodes)
%% Description
% SSATTR single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state
%fixenodes=vector of nodes to keep fixed
%%

stateruns=cell(runs,1);
attrcnt = zeros(1,length(attractors));

global tolerance
global n
tolerance_2 = 10*tolerance ; % should be somewhat higher than attr finding algorithm
newattr=0;
for z = 1:runs

t = 0;
state = initstate;
attr = 0;  
 while (attr == 0)
    t = t + 1;  
    if t > fixedtime
    attr = 1; % you can only reach attractor after perturbation period
    state(1,:) = state(2,:) .* state(3,:); %global activity is multiplication of slow and fast subvariables.
    end
    temp = randperm(n);
    if t <= fixedtime % depending on whether fixed interval passed or not
        for l=1:length(fixednodes)
            temp(temp == fixednodes(l)) = []; % remove the node that should stay fixed
        end
    end
        
    attrf = 0;
% if change = 1 % do fast loop only if there was a change in the slow one
    while (attrf == 0) 
        attrf = 1;
        tmp = randperm(n);
        if t <= fixedtime % depending on whether fixed interval passed or not
            for l=1:length(fixednodes)
                tmp(tmp == fixednodes(l)) = []; % remove the node that should stay fixed
            end
        end
        for j = tmp
            changednode = zupdateas_AC(state,j,'fast'); 
            
        if (abs(state(2,j)-changednode) > tolerance_2)  %&& ~ismember(j,inputindex) % newstate(2,j) ~= changednode;
           state(2,j) = changednode;
           state(1,j) =  state(2,j)*state(3,j); % min(newstate(2,j), newstate(3,j));
           attrf = 0;
%           state(1,inputindex) = input; % keep inputnodes fixed
           break;
        end
        end
    end
%  end
        
    for l = temp
            changednode = zupdateas_AC(state,l,'slow');
         if (abs(state(3,l)-changednode) > tolerance_2) %&& ~ismember(l,inputindex) %newstate(3,l) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
           state(3,l) = changednode;
           state(1,l) = state(2,l)*state(3,l);
           attr = 0;
           break;
         end       
    end

 end
     if ~checkstates(state,attractors,length(attractors))
        disp('new attractor found')
        newattr=newattr+1;
        stateruns{z}=state;
    else
        [~,old] = checkstates(state,attractors,length(attractors));
        attrcnt(old) = attrcnt(old) + 1;
    end  
end
end

function [y2]=stationary(Q)
%Q is the transition matrix of the Markov Chain. We use the method of state
%space reduction to compute the stationary distributions of Q.
    syms z; 
    P=Q;
    [ns ms]=size(P);
    n=ns;
    n8=n;
    while n8>1 
       for n1=1:(n8-1)
              for n2=1:(n8-1)
    P(n1,n2)=P(n1,n2)+P(n1,n8)*P(n8,n2)/(1-P(n8,n8));
              end
       end
    n8=n8-1;
    end

    y=sym(ones(n,1));
       n3=1;
    while n3<n
        x=sym(ones(n3,1));
        for n5=1:n3
        x(n5)=y(n5);
        end        
        z=sum(x.*P(1:n3,n3+1))/(sym(1)-P(n3+1,n3+1));
            y(n3+1)=z/(sym(1)+z);
        for n6=1:n3
            y(n6)=y(n6)*(sym(1)-y(n3+1));
        end
        n3=n3+1;
    end
    y2=double(y);
end

% In silico conditions associated functions:

function [table2save]= Save_insilico_conditions_with_stats(attractors,nodes,inputvalues,initial_attr,sheetNumber,filename,saturation)
% run 3 times the same perturbed conditions in order to compute average and
% standard deviation in state transition.

runs = 3;
results_runs = zeros(runs,length(attractors));
for run = 1 : runs
    run
    [TransitionResult] = inSilicoConditions_return(stauration,attractors,initial_attr,nodes,inputvalues);
    results_runs(run,:) = TransitionResult ;
end

    results_runs(4,:) = mean(results_runs);
    results_runs(5,:) = std(results_runs);

    colNames = {'None','Sox9','Runx2'};
    rowNames = {'run1','run2','run3','average','std'};
    table2save = array2table(results_runs, 'VariableNames', colNames,'RowNames',rowNames);
    
% do not write if ouput was nan (e.g. in case of wrong init_attr indexing):
if ~isnan(mean(results_runs))    
    writetable(table2save,filename,'WriteVariableNames',true,'WriteRowNames',true,'Sheet',sheetNumber,'Range','A2'); 
    xlswrite(filename,[nodes;inputvalues],sheetNumber,'F2');
end

end

function [TransitionResult] = inSilicoConditions_return(saturation,attractors,init_attr,nodes,inputvalues)
%% Description
% Start from the WT attractors and impose perturbed conditions of multiple
% node + study state transition.
%%
    if init_attr> length(attractors)
        TransitionResult = [nan,nan,nan];
        disp('Warning: Wrong indexing of initial state')
        return 
    end
    All_attractors = attractors;
    attr_init = All_attractors{init_attr}; %pick initial state

    runs = 100;
    fixedtime = 800;
    
    global s
    s = zeros(1,5);
    s(1) = 1;
    for i = 2:5
    s(i) = 2* saturation/i;
    end

    per = attr_init;
    per(1,nodes) = inputvalues;% create perturbed state node1
    if ~isequal(attr_init,zupdates_AC(per))
        [attrcnt,~,stateruns] = ssattr(per, runs, All_attractors, fixedtime,nodes);
        TransitionResult = attrcnt;
        storeNewattr=stateruns;
    else
        attrcnt = 0; % all transitions will go to the initial state
        TransitionResult = attrcnt;
        storeNewattr ='no transition';  
    end
                         
end

% pairwise pertrubation screening associated functions:

function [] = Bi_perturbations_noNone(saturation,runumber,attractors,node1,node2,path_dir)
% Bi_perturbations() computes the probability of transition from each
% stable state under single or simultaneous perturbation of 2 nodes
%% input:
%+ staturation factor to balance weight of interactions, default= 2/3
%+ runnumber= the analysis must be run 3 times, indicate the number of
%current run. runnumber value is used to name the saved file with results
%+ attractor_file = name of the .mat output file of the attr finding algorithm
%+ node1 and node2 are the index of nodes to be perturbe

%% output:
% twicepresult: inital attractors in lines and type of perturbation in
% columns. It is a cell array and each cell gives the % of transition to
% other attractors from a given initial state under a given perturbation
% storeNewattr : in some cases, perturbations can lead to a new attractor
% which didn't exist in non perturbed WT system. If it is the case it would be
% reported in this matrix
%%

All_attractors=attractors;
attractors=attractors(2:3); %Only screen Runx2+ and Sox9+ as init state (not robust)

disp (['node' num2str(node1) ' node' num2str(node2) ' started'])

runs = 100;
fixedtime = 1000;

global s
s = zeros(1,5);
s(1) = 1;
for i = 2:5
s(i) = 2* saturation/i;
end

% 4 possible perturbations:
%   (j1,j2)=(0,0)
%   (j1,j2)=(1,1)
%   (j1,j2)=(0,1)
%   (j1,j2)=(1,0)

L=length(attractors);
storeNewattr= cell(L,4); % columns= '(0;0) (1;1) (0;1) (1;0)' 
twicepresult = cell(L,4); % columns= '(0;0) (1;1) (0;1) (1;0)' 

perturb=cell(L,4); %keep track of which node is actually changed in comparison with initial state:
%  - perturb(p,a)= [1,1] if both nodes are changed in perturbation p in comaprison to initial attractor a,
%  - perturb(p,a)=[1,0] or [0,1] if only 1 node changed
%  - perturb(p,a)=[0,0] if no perturbation.

for i = 1:length(attractors)   
    attr = attractors{i}; 
j1 = node1;
j2 = node2;
init=[attr(1,j1),attr(1,j2)];

        for k = 1:2
            per = attr;
            per(1,j1) = k-1;% create perturbed state node1
            per(1,j2) = k-1;% create perturbed state node2
            if ~isequal(attr,zupdates_AC(per))
                [attrcnt,~,stateruns] = attr_transition(per, runs, All_attractors, fixedtime,[j1,j2]);
                twicepresult(i,k) = {attrcnt};
                storeNewattr(i,k)={stateruns};
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                twicepresult(i,k) = {attrcnt};
                storeNewattr(i,k)={'no transition'};  
            end
            
            if init == [per(1,j1),per(1,j2)]
                perturb{i,k}=[0,0]; %no perturbation
            elseif (init(1) == per(1,j1)) && (init(2) ~= per(1,j2))
                perturb{i,k}=[0,1]; %node 1 perturbed only
            elseif (init(1)~=per(1,j1)) && (init(2)==per(1,j2))
                perturb{i,k}=[1,0]; %node 2 perturbed only
            else
                perturb{i,k}=[1,1]; %double perturbation
            end
            %per(1,j1)=k-1; %already set before
            per(1,j2)=~(k-1);
            c=k+2;
            if ~isequal(attr,zupdates_AC(per))
                [attrcnt,~,stateruns] = attr_transition(per, runs, All_attractors, fixedtime,[j1,j2]);
                twicepresult(i,c) = {attrcnt};
                storeNewattr(i,c)={stateruns};  
            else
                attrcnt = zeros(1,length(attractors)); % all transitions will go to the initial state
                twicepresult(i,c) = {attrcnt};
                storeNewattr(i,c)={'no transition'};  
            end 
            
            if init == [per(1,j1),per(1,j2)]
                perturb{i,c}=[0,0]; %no perturbation
            elseif (init(1) == per(1,j1)) && (init(2) ~= per(1,j2))
                perturb{i,c}=[0,1]; %node 1 perturbed only
            elseif (init(1)~=per(1,j1)) && (init(2)==per(1,j2))
                perturb{i,c}=[1,0]; % node 2 perturbed only
            else
                perturb{i,c}=[1,1]; % double perturbation
            end
            
        end
end       
%save variables in file.mat:
string = ['Bi_perturbation' num2str(runumber) '(' num2str(j1) ',' num2str(j2) ')' '.mat'];
f = fullfile(pwd,path_dir,string);
save(f,'saturation','twicepresult','node1','node2','attractors','perturb','storeNewattr')

disp (['node' num2str(node1) ' node' num2str(node2) ' finished'])


end

function [attrcnt,newattr,stateruns] = attr_transition(initstate,runs,attractors, fixedtime, fixednodes)
% single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state
% fixenodes = vector of nodes to keep fixed
% !!  same as SSATTR

stateruns=cell(runs,1);
attrcnt = zeros(1,length(attractors));
global tolerance
global n
tolerance_2 = 10*tolerance; % should be somewhat higher than attr finding algorithm
newattr=0;
for z = 1:runs

t = 0;
state = initstate;
attr = 0;  
 while (attr == 0)
    t = t + 1;
    if t > fixedtime
        attr = 1; % you can only reach attractor after perturbation period
        state(1,:) = state(2,:) .* state(3,:); %global activity is multiplication of slow and fast subvariables.
    end
        temp = randperm(n);
    if t <= fixedtime % depending on whether fixed interval passed or not
        for l=1:length(fixednodes)
            temp(temp == fixednodes(l)) = []; % remove the node that should stay fixed
        end
    end  
    
    attrf = 0;
    while (attrf == 0) 
        attrf = 1;
        tmp = randperm(n);
        if t <= fixedtime % depending on whether fixed interval passed or not
            for l=1:length(fixednodes)
                tmp(tmp == fixednodes(l)) = []; % remove the node that should stay fixed
            end
        end
        for j = tmp
            changednode = zupdateas_AC(state,j,'fast'); 
        if (abs(state(2,j)-changednode) > tolerance_2)
           state(2,j) = changednode;
           state(1,j) =  state(2,j)*state(3,j);
           attrf = 0;
%           state(1,inputindex) = input; % keep inputnodes fixed
           break;
        end
        end
    end
        
    for l = temp
            changednode = zupdateas_AC(state,l,'slow');
         if (abs(state(3,l)-changednode) > tolerance_2)
           state(3,l) = changednode;
           state(1,l) = state(2,l)*state(3,l);
           attr = 0;         
           break;
         end         
    end
 end
     if ~check(state,attractors,length(attractors))
        disp('new attractor found')
        newattr=newattr+1;
        stateruns{z}=state;
    else
        [~,old] = check(state,attractors,length(attractors));
        attrcnt(old) = attrcnt(old) + 1;
     end  
end

function [boolean, number] = check(newstate,state,maxidx)
% function to check whether a state has already been described
boolean = 0;
number = 0;

if maxidx ~= 0
    for i = 1:maxidx
        if  abs(newstate-state{i}) < 10*tolerance 
            boolean = 1;
            number = i;
            
            break;
        end
    end
end
end

end

function [Summary_screening,MyCombi] = screening_summary(C,path_dir)
%% Summarize perturbations screening in big matrix
% output = n-by-a matrix where n is total the number of perturbations a the
% number of attractors. Each value gives the % of transition towards one of
% the target attractors. The initial attractor being the hypertophic one. 

% total number of perturbations= length(combnk(1:60,2))*4 = 7080

%C: full set of possible combinations (e.g. combnk(1:60,2))
nc =size(C,1)*4 ;
 
%%
init_attr=["Sox9","Runx2"];
perturb_type = {'down-down','up-up','down-up','up-down'};
MyCombi=cell(nc,1);
Summary_screening = zeros( nc , 3);
nc=0;

for comb=1:size(C,1) %screen each combination one by one
    node1=C(comb,1);
    node2=C(comb,2);
    string= ['Bi_perturbation1(' num2str(node1) ',' num2str(node2) ').mat'];
    f = fullfile(path_dir,string);
    
    if exist(f,'file') == 2 %if file exists for current combination 
        
        load(f,'twicepresult','attractors','perturb') % load perturbation results of current combination
        
        a = 2;% we are only interested in Runx2+ as starting attractor
       for i=1:4            
           % if perturb ~= [0,0] then perturbation = initial attractor value
           nc = nc+1;
           Summary_screening(nc,:) = twicepresult{a,i};
           
           NameFirst=fromNumToName(node1,'components.mat');
           NameSecond=fromNumToName(node2,'components.mat');
           MyCombi{nc} = [char(NameFirst) '-' char(NameSecond) ' ' perturb_type{i}];     
           %MyCombi{nc} = [char(NameFirst) ' - ' char(NameSecond) ' ' num2str(i)]  
       end
             
    else   
       disp(['combination ' num2str(C(comb,1)) ',' num2str(C(comb,2)) ' doesn''t exist!'])
    end
end

save('Summary_screening_from_Hyp.mat','Summary_screening','MyCombi')
end

function[VarName]=fromNumToName(VarNum,ref)
%% Description
%This function translates variable indexes into their biological names
% Varnum : int array
% VarName: string array
%ref: string giving the name of the matrice to load (ex: 'componentnames.mat'
%%

%load the list of component in which idx corresponds to the variable number
load(ref,'componentnames'); %structure with string array.

componentnames=string(componentnames(1,:));
ntot=length(componentnames);
sz=size(VarNum);

%initialisation of the output string array:
VarName=strings(sz);

%Browse the input variable values:
for i=1:sz(1)
    for j=1:sz(2)
        %make sure the input variable is within the scope of number of
        %existing component inmodel (ntot=60 components)
        if (VarNum(i,j)>0 && VarNum(i,j)<=ntot) 
            %get the corresponding variable name
            VarName(i,j)=componentnames(VarNum(i,j)); 
        elseif (VarNum(i,j)>ntot+1 && VarNum(i,j)<=ntot*2) 
            %get the corresponding variable name
            VarName(i,j)=componentnames(VarNum(i,j)-46); 
        else
            VarName(i,j)=0; % set 0 if the input is not valid
        end
    end
end

end

function [h]=plot_PCA_pertub_screening(Summary_screening,MyCombi)
%% Principal Component Analysis of perturbation screening
%PCA on transitions percentage results and plot it

%Summary_screening : each column is an attractor, each row a perturbation.
%values give the % of transition towards the respective attractors. 
%%

% PCA with attractors as variables and  perturbations as observation:
[coeff,score,latent,~,explained] = pca(Summary_screening);

%%
%Colored biplot per threshold, 3 groups:
    % < 70% Healthy
    % >= 70% Healthy
    % >=70% healthy but <5% None

% Store handle to biplot
    h =biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'None','Healthy','Hyp'});

% Identify each handle
    hID = get(h, 'tag');

% Isolate handles to scatter points
    hPt = h(strcmp(hID,'obsmarker'));

% Identify groups according to threshold
    threshold = 70;
    indx2 = find(Summary_screening(:,2)>= threshold);
    indx3 = find(Summary_screening(:,2)>= threshold & Summary_screening(:,1)<= 5);
    grp = ones(length(Summary_screening),1);
    grp(indx2)=2;
    grp(indx3)=3;
    grpID ={'<70% Healthy','>= 70% Healthy','>= 70% Healthy & <5% None'}

% Assign colors and legend display name
%clrMap = lines(length(unique(grp)));   % using 'lines' colormap
    clrMap = [0.6,0.6,0.6 ; 0,0.4470,0.7410 ; 0.8500,0.3250,0.0980];
    for i = 1:max(grp)
        set(hPt(grp==i), 'Color', clrMap(i,:), 'DisplayName', grpID{i})
    end

% add legend to identify cluster
    [~, unqIdx] = unique(grp);
    legend(hPt(unqIdx))


xlabel(sprintf('Component 1 (%0.1f %%)',explained(1)))
ylabel(sprintf('Component 2 (%0.1f %%)',explained(2)))
end


function [MyCombi]=select_perturbations(C,init_attr,percentSox9,percentRunx2,path_dir)
%Select the combinations the most relevant for in vitro testing
MyCombi=struct;
j=0;
    for comb=1:size(C,1) %screen each combination one by one
        node1=C(comb,1);
        node2=C(comb,2);
        string= ['Bi_perturbation1(' num2str(node1) ',' num2str(node2) ').mat'];
        f = fullfile(path_dir,string);
        if exist(f,'file') == 2 %if file exists for current combination 
            twicepresult={};
            attractors={};
            load(f,'twicepresult','attractors') % load perturbation results of current combination
           for a=1:length(init_attr) %from Sox9 (1) and from Runx2 (2)
    
            perturbation = ["down-down","up-up","down-up","up-down"];
            for i=1:4  
                percent=twicepresult{a,i}; %twicepresult{f,i};
                NameFirst=fromNumToName(node1,'components.mat');
                NameSecond=fromNumToName(node2,'components.mat');

                if (init_attr(a)=="Sox9" && percentSox9 <= 50)            
                    if (percent(2)<=percentSox9 && percent(3)>=percentRunx2)
                        j=j+1;
                        %at least x% of perturbations lead to Runx2
                        MyCombi(j).Combi=[node1 node2]; %store the nodes
                        %1: down-down, 2: up-up, 3:down-up, 4: up-down:
                          
                        MyCombi(j).Perturb= perturbation(i); % store the type of perturbation
                        MyCombi(j).init_attractor=init_attr(a);
                        MyCombi(j).initvalues=[attractors{a}(1,node1) attractors{a}(1,node2)]; %valie of node in initila attractor
                        MyCombi(j).Names=[char(NameFirst) ' - ' char(NameSecond)];

                    end               
                elseif (init_attr(a)=="Runx2" && percentRunx2 <= 50)
                    if (percent(2)>=percentSox9 && percent(3)<=percentRunx2)
                        j=j+1;
                        %at least 50% of perturbations lead to Sox9
                        MyCombi(j).Combi=[node1 node2]; %store the nodes
                        MyCombi(j).Perturb=perturbation(i); % store the type of perturbation %1 for down-down, 2 for up-up, 3 for down-up, 4 for down-up
                        MyCombi(j).init_attractor=init_attr(a);
                        MyCombi(j).initvalues=[attractors{a}(1,node1) attractors{a}(1,node2)]; %valie of node in initila attractor
                        MyCombi(j).Names=[char(NameFirst) ' - ' char(NameSecond)];
                        
                    end
                end
            end
           end

        else
            disp(['combination ' num2str(C(comb,1)) ',' num2str(C(comb,2)) ' doesn''t exist!'])
        end  
    end
end

function [report]=Report_selected_conditions()
report = struct();

load('Combi2test100_0.mat','SelectedCombi') 
combi_100=struct2table(SelectedCombi); % 100% transition
writetable(combi_100,'CombinationSelection.xls')
report.full_transition = combi_100;

load('Combi2test90_10.mat','SelectedCombi')
combi_90=struct2table(SelectedCombi); % >90% transition

load('Combi2test80_20.mat','SelectedCombi')
combi_80=struct2table(SelectedCombi); % >80% transition

load('Combi2test70_30.mat','SelectedCombi')
combi_70=struct2table(SelectedCombi); % >70% transition



%Find conditions in 90% threshold that are not in 100% threshold, etc..
    result1 = setdiff(combi_90,combi_100);
    writetable(result1,'CombinationSelection.xls','Sheet','Sheet2')

    result2 = setdiff(combi_80,combi_90);
    writetable(result2,'CombinationSelection.xls','Sheet','Sheet3')

    result3 = setdiff(combi_70,combi_80);
    writetable(result3,'CombinationSelection.xls','Sheet','Sheet4')

report.moreThan90 = result1; %> 90% but <100%
report.moreThan80 = result2; %> 80% but <90%
report.moreThan70 = result3; %> 70% but <80%

end
