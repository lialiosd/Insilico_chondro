function [TransitionResult] = inSilicoConditions_return(saturation,attractor_file,init_attr,nodes,inputvalues)
%% Description
% Start from the WT attractors and impose perturbed conditions of multiple
% node + study state transition.
%%
attractorsload = load(attractor_file,'attractors'); 
    All_attractors = attractorsload.attractors;
    attr_init = All_attractors{init_attr}; %pick initial state

    runs = 100;
    fixedtime = 800;
    
    global n
    n = 60;   %number of nodes in current system
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

function [attrcnt,newattr,stateruns] = ssattr(initstate,runs,attractors, fixedtime, fixednodes)
%% Description
% SSATTR single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state
%fixenodes=vector of node to keep fixed
%%

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
            
        if (abs(state(2,j)-changednode) > tolerance)  %&& ~ismember(j,inputindex) % newstate(2,j) ~= changednode;
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
         if (abs(state(3,l)-changednode) > tolerance) %&& ~ismember(l,inputindex) %newstate(3,l) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
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
        if  abs(newstate-state{i}) < 10*tolerance % isequal(newstate, state{i});
            boolean = 1;
            number = i;
            
            break;
        end
    end
end
end
