%% Summarize perturbations screening in big matrix
% output = n-by-a matrix where n is total the number of perturbations a the
% number of attractors. Each value gives the % of transition towards one of
% the target attractors. The initial attractor being the hypertophic one. 

% total number of perturbations= length(combnk(1:60,2))*4 = 7080

C=combnk(1:60,2); %create the full set of possible combinations
n =length(combnk(1:60,2))*4 ;
 
%%
init_attr=["Sox9","Runx2"];
perturb_type = {'down-down','up-up','down-up','up-down'};
MyCombi=cell(n,1);
Summary_screening = zeros( n , 3);
n=0;

for comb=1:length(C) %screen each combination one by one
    node1=C(comb,1);
    node2=C(comb,2);
    string= ['Bi_perturbation1(' num2str(node1) ',' num2str(node2) ').mat'];

    if exist(string,'file') == 2 %if file exists for current combination 
        
        load(string,'twicepresult','attractors','perturb') % load perturbation results of current combination
        
        a = 2;% we are only interested in Runx2+ as starting attractor
       for i=1:4            
           % if perturb ~= [0,0] then perturbation = initial attractor value
           n = n+1
           Summary_screening(n,:) = twicepresult{a,i};
           
           NameFirst=fromNumToName(node1,'component.mat');
           NameSecond=fromNumToName(node2,'component.mat');
           MyCombi{n} = [char(NameFirst) '-' char(NameSecond) ' ' perturb_type{i}];     
           %MyCombi{n} = [char(NameFirst) ' - ' char(NameSecond) ' ' num2str(i)]  
       end
             
    else   
       disp(['combination ' num2str(C(comb,1)) ',' num2str(C(comb,2)) ' doesn''t exist!'])
    end
end

save('Summary_screening_from_Hyp.mat','Summary_screening','MyCombi')

%% Required functions
function[VarName]=fromNumToName(VarNum,ref)
%% Description
%This function translate variable numbers into their biological names
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