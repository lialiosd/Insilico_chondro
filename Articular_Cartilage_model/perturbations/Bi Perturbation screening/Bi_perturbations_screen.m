%Perform perturbation of 2 nodes at the same time, screening of all
%possible combination run in parallel (for HPC) 

C=combnk(1:60,2); % all combinations of 2 among 60. 

disp(['matlab number of cores : '   num2str(feature('numCores'))])

parfor comb=1:length(C)
    Ccomb=C(comb,:)
    nod1=Ccomb(1)
    nod2= Ccomb(2)
	disp(comb)
        Bi_perturbations_noNone(2/3,1,'attractorAC_WT.mat',nod1,nod2)    
end

function [] = Bi_perturbations_noNone(saturation,runumber,attractor_file,node1,node2)
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
attractorsload=load(attractor_file,'attractors'); 
attractors=attractorsload.attractors;
All_attractors=attractors;
attractors=attractors(2:3); %Only screen Runx2+ and Sox9+ as init state (not robust)

disp (['node' num2str(node1) ' node' num2str(node2) ' started'])
tic

runs = 100;
fixedtime = 1000;
global n
n = 60;   %number of nodes in current system
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
                [attrcnt,~,stateruns] = ssattr(per, runs, All_attractors, fixedtime,[j1,j2]);
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
                [attrcnt,~,stateruns] = ssattr(per, runs, All_attractors, fixedtime,[j1,j2]);
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
save(string,'saturation','twicepresult','node1','node2','attractors','perturb','storeNewattr')
   
disp (['node' num2str(node1) ' node' num2str(node2) ' finished'])
toc


end

function [attrcnt,newattr,stateruns] = ssattr(initstate,runs,attractors, fixedtime, fixednodes)
% SSATTR single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state
% fixenodes = vector of node to keep fixed

stateruns=cell(runs,1);
attrcnt = zeros(1,length(attractors));
global tolerance
global n
tolerance = 1e-2; % should be somewhat higher than attr finding algorithm
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
        if (abs(state(2,j)-changednode) > tolerance)
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
         if (abs(state(3,l)-changednode) > tolerance)
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

function [boolean, number] = checkstates(newstate,state,maxidx)
% function to check whether a state has already been described
global tolerance
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
