%% Principal Component Analysis of perturbation screening
%Description : load petrubation screening summary and PCA on transitions
%percentage results. 
%%
load('Summary_screening_from_Hyp.mat','Summary_screening','MyCombi')
%Summary_screening : each column is an attractor, each row a perturbation.
%values give the % of transition towards the respective attractors. 

% PCA with attractors as variables and  perturbations as observation
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