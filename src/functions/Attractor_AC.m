function [d, sizes, da, percent, attractors] = Attractor_AC(nstates, saturation,input, inputindex)
%% Algorithm to compute the model stable states through a Monte Carlo method
% This file is part of the Insilico_chondo repository
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2021 - KU Leuven

%File author(s): Raphaëlle Lesage (partially based on J.Kerkhofs et al. PLoS
%One (2016))

%Distributed under the GPLv3 License.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>

%%
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

% rng('shuffle')
global n
n = 60;
global tolerance
tolerance = 1e-3;

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
    ptr
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
string = 'attractorAC_WT'; 

[attractors_read,T_sizes,T_percent,Attr_profiles] = reader_friendly_output(sizes, percent, attractors);

save([string,'.mat'], 'd', 'sizes', 'da', 'percent', 'attractors','attractors_read');
save([string,'_readerFriendly.mat'], 'T_sizes','T_percent','Attr_profiles');

end

%% Functions

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