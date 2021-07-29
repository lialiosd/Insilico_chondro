function [perresult, stabattr, nb_transition,storeNewattr,transition,twicepresult,distribution] = StateStab_perturbation(saturation,runnumber,attractor_file)
% StateStab_perturbation: Calculates the stability of a states by perturbing network nodes one by
% one

[attractors] = format_attractors(attractor_file);
disp(['formatting done...'])

%set parameters:
runs = 100; %number of repetition of the perturbation from init-state
fixedtime = 1000; %duration of perturbation (in time steps)

if ~isnumeric(saturation)
   saturation = str2num(saturation);
end

global n
n = 60;   
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

% save variables in file.mat:
string = ['SingleNodePertub' num2str(runnumber) '_sat' num2str(saturation) '.mat'];
save(string, 'stabattr', 'perresult', 'trlik', 'trlik2','saturation','storeNewattr','transition','distribution')
   

end

function [attrcnt,newattr,stateruns] = ssattr(initstate,runs,attractors, fixedtime, fixednode)
% SSATTR single state attractor analysis: investigate probability of ending
% up in a certain attractor from a given starting state

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
    state(1,:) = state(2,:) .* state(3,:);
    end
    temp = randperm(n);
    if t <= fixedtime % depending on whether fixed interval passed or not
    temp(temp == fixednode) = []; % remove the node that should stay fixed
    end
        
    attrf = 0;
% if change = 1 % do fast loop only if there was a change in the slow one
    while (attrf == 0) 
        attrf = 1;
        tmp = randperm(n);
        if t <= fixedtime % depending on whether fixed interval passed or not
        tmp(tmp == fixednode) = []; % remove the node that should stay fixed
        end
        for j = tmp
            changednode = zupdateas_AC(state,j,'fast'); 
        if (abs(state(2,j)-changednode) > tolerance)  %&& ~ismember(j,inputindex) % newstate(2,j) ~= changednode; % 1st row contains the effective values, the second contains the protein values, the third mRNA
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
%            change = 1;
           
            %          state(1,inputindex) = input; % keep inputnodes fixed
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
        if  abs(newstate-state{i}) < 100*tolerance % isequal(newstate, state{i});
            boolean = 1;
            number = i;
            
            break;
        end
    end
end
end

function [attractors] = format_attractors(attractor_file)
% Transform cell with node names into matrix with numerical values only
    attractorsload = load(attractor_file,'attractors_read'); 
    attractors = attractorsload.attractors_read;
    %% if nodenames
    for i=1:length(attractors)
        attractors{i}=attractors{i}(2:4,:);
    end
    attractors = cellfun(@str2double,attractors,'UniformOutput',false);
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