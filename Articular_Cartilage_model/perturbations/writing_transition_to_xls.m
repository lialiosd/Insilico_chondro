
%% Write the result of 'perresult' in xls sheet to analyse the transition per node
%NB: for run2 and run3: change the sheet number to avoid overwritting

%% License
% This file is part of the Insilico_chondo repository
% (https://github.com/Rapha-L) 

%Copyright (c) 2017-2021 - KU Leuven

%File author(s): Raphaëlle Lesage (contact: liesbet.geris@kuleuven.be)

%Distributed under the GPLv3 License.
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>


%% definition of variables and parameters
n=60; %number of nodes
a=2; % number of the starting attractor for the transition (e.g. SOX9+=attr2)

%% loading of file
%load the mat containing 'perresult' :
load('SingleNodePertub1sat0.66667.mat') %e.g. 'SingleNodePertub.....mat'

%load variable names:
load('components.mat','componentnames');
names=[componentnames{1,:}];

%name of the excel file? where data are written and saved
xl_file= 'transition_pernode_AC_fromSox9.xls';


%% wrinting in excel

for c=[1,2]
    if (c==1)
        perturbation_type="d";
    elseif (c==2)
        perturbation_type="u";
    end
    sheet=c; % c for run1 %c+2 for run2  %c+4 for run3
    xlswrite(xl_file,perturbation_type,sheet,'A1')
    xlswrite(xl_file,["attr1", "attr2", "attr3"],sheet,'B2')
    
    %write a colonne with the name of node in excel sheet
    xlswrite(xl_file,names',sheet,'A3');
    for i=1:n
        j=i+2; % keep few lines for headers above the table
        
        %write value of transition from attractor 'a' to others for each node
        %for perturbation 'c'
        if (~isempty(perresult{a,i,c}))
            xlswrite(xl_file,perresult{a,i,c},sheet,['B' num2str(j)]);
        else
             xlswrite(xl_file,[0 0 0],sheet,['B' num2str(j)]);
        end
        
    end
end