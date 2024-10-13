
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILES UTILIZED: FinalData (folder), Results_DMM-Mix (folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Variable Workspace and Command Window
clear
clc
close all force

mkdir       Results_DMM-Mix
outDir      = 'Results_DMM-Mix'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Parameters

vThresh     = 1;  % Variance Threshold for Mouse PCs
qThresh     = 0.2; % Mouse Gene Selction q Threshold   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ DATA IN

dmm_dataIn      = readtable('FinalData/GSE41342_DMM_mouse_homologGenes.txt','Delimiter','\t','ReadRowNames',1);
dmm_pheno       = readtable('FinalData/GSE41342_DMM_mouse_phenotypes.txt','Delimiter','\t','ReadRowNames',0);

% drop non-Sham controls
dmm_dataIn(:,[1:3])     = [];
dmm_pheno([1:3],:)      = [];

oa999_dataIn    = readtable('FinalData/GSE117999_human_genes.txt','Delimiter','\t','ReadRowNames',1);
oa999_pheno     = readtable('FinalData/GSE117999_human_phenotypes.txt','Delimiter','\t','ReadRowNames',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DMM Mouse Processing

mouse_In    = dmm_dataIn;
mouse_phen  = dmm_pheno;

% Create Mouse DEG Testing Tables
T                                   = splitvars(table(mouse_phen.Surgery_DMM_is_1,...
                                                mouse_phen.Time,zscore(mouse_In{:,:}')));
T.Properties.VariableNames(1:2)     = {'Surgery'; 'Time'};
T.Properties.VariableNames(3:end)   = mouse_In.Properties.RowNames;

% Get model statistics table dimensions
mdl_tmp                             = fitglm(T(:,1:3));
pArray                              = zeros(height(mouse_In),1+height(mdl_tmp.Coefficients));

% Light feature selection, variable mouse genes 
    parfor i = 1:height(mouse_In)
        mdl             = fitglm(T(:,[1 2 i+2]));
        pArray(i,:)     = [mdl.coefTest mdl.Coefficients.pValue'];
    end
    
% Output statistics
pTable                                  = splitvars(table(pArray,'RowNames',mouse_In.Properties.RowNames));
pTable.Properties.VariableNames(1)      = {'model'};
pTable.Properties.VariableNames(2:end)  = mdl_tmp.Coefficients.Properties.RowNames;
[~,k]                                   = size(pTable);
for i = 1:k
    qArray(:,i) = mafdr(pArray(:,i),'BHFDR','True');
end
qTable      = pTable;
qTable{:,:} = qArray;

% Output Feature Selection Statistics
writetable(pTable,[outDir '/DMM_featureSelection_pValues.txt'],'WriteRowNames',1)
writetable(qTable,[outDir '/DMM_featureSelection_qValues.txt'],'WriteRowNames',1)

% Pick DMM Injury/Disease-associated Features, Filter Mouse/Human Datasets for genes
idx                 = find(qTable{:,3} < qThresh);
mouse_data          = mouse_In(idx,:);
[~,h1,h2]           = intersect(oa999_dataIn.Properties.RowNames,oa999_dataIn.Properties.RowNames);
oa999_data          = oa999_dataIn(h2,:);
[~,ih,im]           = intersect(oa999_data.Properties.RowNames,mouse_data.Properties.RowNames);
oa999_data          = oa999_data(ih,:);
mouse_data          = mouse_data(im,:);

if ~isequal(oa999_data.Properties.RowNames,mouse_data.Properties.RowNames);
    disp('!!ERROR!! DMM - OA999 Mismatch');end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DMM Mouse PCA and PC-R Model Building and Analysis

% Mouse PCA model for making the scores plots
[~, ~, ~, ~, EXP]   = pca(zscore(mouse_data{:,:}'));
numPC               = min(find(EXP < vThresh)) - 1; % Variance Threshold
[XL, XS, ~, ~, EXP] = pca(zscore(mouse_data{:,:}'),'NumComponents',numPC);

%Save loadings with gene names for GSEA
XL_genes = array2table(XL);
XL_genes.Properties.RowNames = mouse_data.uniqueIDs;
filename = 'DMMMouseLoadings_w_genes.xlsx';
writetable(XL_genes, filename, 'Sheet', 'Sheet1','WriteRowNames', true);

% Mouse Data PC-R Model
mouse_mdlTable       = splitvars(table(XS,mouse_phen.Time,mouse_phen.Surgery_DMM_is_1));
mouse_mdlTable.Properties.VariableNames(numPC+1:end) = {'Time';'Disease'};

% Train PC-R, Mouse PC's predict Mouse Conditions/Phenotype
mdl_mouseTime        = fitglm(mouse_mdlTable(:,[1:end-1]));
mdl_mouseOA          = fitglm(mouse_mdlTable(:,[1:end-2 end]));

% Incorporating the Time covariate 
mdl_specMouse = 'Disease ~ XS_1 + XS_2 + XS_6 + XS_7 + XS_1*Time + XS_2*Time + XS_6*Time + XS_7*Time';
mdl_mouseFull = fitglm(mouse_mdlTable,mdl_specMouse);

% Output DMM Loadings and PC-R Statistics For Downstream Analysis
loadingTable                            = splitvars(table(XL,'RowNames',mouse_data.Properties.RowNames));
writetable(loadingTable ,[outDir '/DMM_mouseGeneLoadings.txt'],'WriteRowNames',1)

writetable(mdl_mouseTime.Coefficients,[outDir '/DMM_mouse_TimePC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/DMM_mouse_TimePC-R_pValue.csv'],mdl_mouseTime.coefTest)

writetable(mdl_mouseOA.Coefficients,[outDir '/DMM_mouse_DiseasePC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/DMM_mouse_DiseasePC-R_pValue.csv'],mdl_mouseOA.coefTest)

writetable(mdl_mouseFull.Coefficients,[outDir '/DMM_mouse_Disease_x_Time_PC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/DMM_mouse_Disease_x_Time_PC-R_pValue.csv'],mdl_mouseFull.coefTest)

% Percent Variance Explained By Mouse PCs
for i = 1:numPC
    exp_oa999(i,1)      =   XL(:,i)'*zscore(oa999_data{:,:}')'*zscore(oa999_data{:,:}')*XL(:,i)./...
        sum(diag(XL'*zscore(oa999_data{:,:}')'*zscore(oa999_data{:,:}')*XL));
end

% Variance Explained In Mouse and Human Model
vTable = splitvars(table([EXP(1:numPC), 100*exp_oa999]));
vTable.Properties.VariableNames = {'Mouse';'OA999'};
disp(vTable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merge Human Data Projections

human_data1      = oa999_data;
human_phen1      = oa999_pheno;

pc              = [zscore(human_data1{:,:}')*XL];  %human data multiplied by mouse loadings
age             = [human_phen1.AgeNumeric];
sex             = [human_phen1.SexNumeric];
oa              = [human_phen1.Disease];
batch           = [-1*ones(height(human_phen1),1)];

h_mdlTable                                          = splitvars(table(pc,age,sex,batch,oa));
h_mdlTable.Properties.VariableNames(numPC+1:end)    = {'Age';'Sex';'Batch';'Disease'};

% Human Models (GLM) used to pick PCs for the TransComp-R model
mdl_OA         = fitglm(h_mdlTable(:,[1:numPC end]));
mdl_sex        = fitglm(h_mdlTable(:,[1:numPC end-2]));
mdl_age        = fitglm(h_mdlTable(:,[1:numPC end-3]));


% Human TransComp-R Models (Combining Mouse PC's with Human Covariates) N = 58 human samples
pc_spec     = 'Disease ~ pc_4 + pc_8 + pc_10 + pc_13';
main_spec   = 'Disease ~ pc_4 + pc_8 + pc_10 + pc_13 + Age + Sex';
age_spec    = 'Disease ~ pc_4 + pc_8 + pc_10 + pc_13 + Age + Sex + pc_4*Age + pc_8*Age + pc_10*Age + pc_13*Age';
sex_spec    = 'Disease ~ pc_4 + pc_8 + pc_10 + pc_13 + Age + Sex + pc_4*Sex + pc_8*Sex + pc_10*Sex + pc_13*Sex';
full_spec   = 'Disease ~ pc_4 + pc_8 + pc_10 + pc_13 + Age + Sex + pc_4*Age + pc_8*Age + pc_10*Age + pc_13*Age + pc_4*Sex + pc_8*Sex + pc_10*Sex + pc_13*Sex';

c_mdl_pcs       = fitglm(h_mdlTable,pc_spec);
c_mdl_main      = fitglm(h_mdlTable,main_spec);
c_mdl_age       = fitglm(h_mdlTable,age_spec);
c_mdl_sex       = fitglm(h_mdlTable,sex_spec);
c_mdl_full      = fitglm(h_mdlTable,full_spec);

writetable(c_mdl_pcs.Coefficients,[outDir '/mergedHumans_mousePCOnly_TransComp-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/mergedHumans_mousePCOnly_TransComp-R_pValue.csv'],c_mdl_pcs.coefTest)
writetable(c_mdl_full.Coefficients,[outDir '/mergedHumans_full_TransComp-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/mergedHumans_full_TransComp-R_pValue.csv'],c_mdl_full.coefTest)

slim_spec     = 'Disease ~ pc_2 + Age + Sex + pc_2*Age + pc_2*Sex';
c_mdl_slim    = fitglm(h_mdlTable,slim_spec);


%% Figure Production

% Create groups of indices for various phenotypes
[APM]     = find(oa999_pheno.Disease == -1);
[OA]     = find(oa999_pheno.Disease == 1);
[dmm]      = find(mouse_mdlTable.Disease == 1);
[sham]    = find(mouse_mdlTable.Disease == -1);

% FIGURE 5B
% Plot human PC scores in mouse space : PC8 vs. PC13 colored by disease
figure(1)
hold on
scatter( pc(APM, 8), pc(APM, 13),135, 'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75);
scatter( pc(OA, 8), pc(OA, 13), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75); 
xlabel('PC8') ; 
ylabel('PC13') ; 
title('DMM OA vs. APM Human Scores in Mouse Space') ; 
legend('APM', 'OA');
xlim([-4 3.5]); 
hold off

% FIGURE 5C
% Plot human PC scores in mouse space : Human Age vs. PC8 colored by disease
figure(2)
hold on 
scatter( age(APM,:), pc(APM, 8),135, 'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75);
scatter( age(OA,:), pc(OA, 8), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75);
xlabel('Age (years)') ; 
ylabel('PC8') ; 
title('DMM OA vs. APM Human Scores in Mouse Space') ; 
legend('APM', 'OA');
xlim([25 85]); 
hold off

% FIGURE 3C
% Plot mouse PC1 vs. PC2 color coded by mouse injury status 
figure(3)
hold on
scatter(mouse_mdlTable.XS_1(sham), mouse_mdlTable.XS_2(sham), 135, 'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75) ; 
scatter(mouse_mdlTable.XS_1(dmm), mouse_mdlTable.XS_2(dmm), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75); 
xlabel('XS1') ; 
ylabel('XS2') ; 
title('DMM Mouse Linear Regression PCs by Injury ') ; 
legend('Sham', 'DMM');
xlim([-50 70]);  
axis equal;
hold off
