%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILES UTILIZED: FinalData (folder), Results_ACL-Mix (folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Variable Workspace and Command Window
clear
clc
close all force

mkdir       Results_ACL-Mix
outDir      = 'Results_ACL-Mix'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Parameters

vThresh     = 1;  % Variance Threshold mouse PCs
qThresh     = 0.2; % mouse gene selction q Threshold   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% READ DATA IN

acl_dataIn      = readtable('FinalData/GSE112641_mouse_homologGenes.txt','Delimiter','\t','ReadRowNames',1);
acl_pheno       = readtable('FinalData/GSE112641_mouse_phenotypes.txt','Delimiter','\t','ReadRowNames',0);

% drop non BL/6 mice 
keep_idx = find(acl_pheno.StrainNumeric == 2);
acl_dataIn = acl_dataIn(:,keep_idx);
acl_pheno = acl_pheno(keep_idx,:);

oa999_dataIn    = readtable('FinalData/GSE117999_human_genes.txt','Delimiter','\t','ReadRowNames',1);
oa999_pheno     = readtable('FinalData/GSE117999_human_phenotypes.txt','Delimiter','\t','ReadRowNames',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACLR Mouse Processing

mouse_In    = acl_dataIn;
mouse_phen  = acl_pheno;

% Create Mouse DEG Testing Tables
T                                   = splitvars(table(mouse_phen.InjuryNumeric,...
                                                mouse_phen.SampleCollectionNumeric_Days,zscore(mouse_In{:,:}')));
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
writetable(pTable,[outDir '/ACL_featureSelection_pValues.txt'],'WriteRowNames',1)
writetable(qTable,[outDir '/ACL_featureSelection_qValues.txt'],'WriteRowNames',1)

% Pick ACL Injury/Disease-associated Features, Filter Mouse/Human Datasets for genes
idx                 = find(qTable{:,3} < qThresh);
mouse_data          = mouse_In(idx,:);
[~,h1,h2]           = intersect(oa999_dataIn.Properties.RowNames,oa999_dataIn.Properties.RowNames);
oa999_data          = oa999_dataIn(h2,:);

if ~isequal(oa999_data.Properties.RowNames,oa999_data.Properties.RowNames);
    disp('!!ERROR!! Human Genes Mismatch');end

[~,ih,im]           = intersect(oa999_data.Properties.RowNames,mouse_data.Properties.RowNames);
oa999_data          = oa999_data(ih,:);
mouse_data          = mouse_data(im,:);

if ~isequal(oa999_data.Properties.RowNames,mouse_data.Properties.RowNames);
    disp('!!ERROR!! ACL - OA999 Mismatch');end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACLR Mouse PCA and PC-R Model Building and Analysis
 
% Mouse PCA model for making the scores plots 
[~, ~, ~, ~, EXP]   = pca(zscore(mouse_data{:,:}')); 
numPC               = min(find(EXP < vThresh)) - 1; % Variance Threshold
[XL, XS, ~, ~, EXP] = pca(zscore(mouse_data{:,:}'),'NumComponents',numPC);

%Save loadings with gene names for GSEA
XL_genes = array2table(XL);
XL_genes.Properties.RowNames = mouse_data.uniqueIDs;
filename = 'ACLRMouseLoadings_w_genes.xlsx';
writetable(XL_genes, filename, 'Sheet', 'Sheet1','WriteRowNames', true);

% Mouse Data PC-R Model
mouse_mdlTable       = splitvars(table(XS,mouse_phen.SampleCollectionNumeric_Days,mouse_phen.InjuryNumeric));
mouse_mdlTable.Properties.VariableNames(numPC+1:end) = {'Time';'Disease'};

% Train PC-R, Mouse PCs predict Mouse Conditions/Phenotype
mdl_mouseTime        = fitglm(mouse_mdlTable(:,[1:end-1]));
mdl_mouseOA          = fitglm(mouse_mdlTable(:,[1:end-2 end]));

% Incorporating the Time covariate
mdl_specMouse = 'Disease ~ XS_1 + XS_2 + XS_3 + XS_4 + XS_5 + XS_6 + XS_7 + XS_1*Time + XS_2*Time + XS_3*Time + XS_4*Time + XS_5*Time + XS_6*Time + XS_7*Time';
mdl_mouseFull = fitglm(mouse_mdlTable,mdl_specMouse);

% Output ACL Loadings and PC-R Statistics For Downstream Analysis
loadingTable                            = splitvars(table(XL,'RowNames',mouse_data.Properties.RowNames));
writetable(loadingTable ,[outDir '/ACL_mouseGeneLoadings.txt'],'WriteRowNames',1)

writetable(mdl_mouseTime.Coefficients,[outDir '/ACL_mouse_TimePC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/ACL_mouse_TimePC-R_pValue.csv'],mdl_mouseTime.coefTest)

writetable(mdl_mouseOA.Coefficients,[outDir '/ACL_mouse_DiseasePC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/ACL_mouse_DiseasePC-R_pValue.csv'],mdl_mouseOA.coefTest)

writetable(mdl_mouseFull.Coefficients,[outDir '/ACL_mouse_Disease_x_Time_PC-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/ACL_mouse_Disease_x_Time_PC-R_pValue.csv'],mdl_mouseFull.coefTest)

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

pc              = [zscore(human_data1{:,:}')*XL];
age             = [human_phen1.AgeNumeric];
sex             = [human_phen1.SexNumeric];
oa              = [human_phen1.Disease];
batch           = [-1*ones(height(human_phen1),1)];

h_mdlTable                                          = splitvars(table(pc,age,sex,batch,oa));
h_mdlTable.Properties.VariableNames(numPC+1:end)    = {'Age';'Sex';'Batch';'Disease'};

% Human Models (GLM) used to pick PCs for the TransComp-R model
mdl_OA         = fitglm(h_mdlTable(:,[1:numPC end]));
mdl_set        = fitglm(h_mdlTable(:,[1:numPC end-1]));
mdl_sex        = fitglm(h_mdlTable(:,[1:numPC end-2]));
mdl_age        = fitglm(h_mdlTable(:,[1:numPC end-3]));

% Human TransComp-R Models (Combining Mouse PC's with Human Covariates) N = 58 human samples

pc_spec     = 'Disease ~ pc_3 + pc_4 + pc_6';
main_spec   = 'Disease ~ pc_3 + pc_4 + pc_6 + Age + Sex';
age_spec    = 'Disease ~ pc_3 + pc_4 + pc_6 + Age + Sex + pc_3*Age + pc_4*Age + pc_6*Age';
sex_spec    = 'Disease ~ pc_3 + pc_4 + pc_6 + Age + Sex + pc_3*Sex + pc_4*Sex + pc_6*Sex';
full_spec   = 'Disease ~ pc_3 + pc_4 + pc_6 + Age + Sex + pc_3*Age + pc_4*Age + pc_6*Age + pc_3*Sex + pc_4*Sex + pc_6*Sex';
simpleSpec = 'Disease ~ pc_3 + pc_4 + pc_6 + Age + Sex + pc_3*Age + pc_4*Age + pc_6*Age + pc_3*Sex + pc_4*Sex + pc_6*Sex';

% Regression on most complex model
c_mdl_simple = fitglm(h_mdlTable,simpleSpec) ; 

% Regression model using each specification 
c_mdl_pcs       = fitglm(h_mdlTable,pc_spec);
c_mdl_main      = fitglm(h_mdlTable,main_spec); 
c_mdl_age       = fitglm(h_mdlTable,age_spec);
c_mdl_sex       = fitglm(h_mdlTable,sex_spec);
c_mdl_full      = fitglm(h_mdlTable,full_spec); 

writetable(c_mdl_pcs.Coefficients,[outDir '/mergedHumans_mousePCOnly_TransComp-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/mergedHumans_mousePCOnly_TransComp-R_pValue.csv'],c_mdl_pcs.coefTest)
writetable(c_mdl_full.Coefficients,[outDir '/mergedHumans_full_TransComp-R.txt'],'WriteRowNames',1)
csvwrite([outDir '/mergedHumans_full_TransComp-R_pValue.csv'],c_mdl_full.coefTest)

%% Figure Production

% Create groups of indices for various phenotypes 
[inj]      = find(acl_pheno.InjuryNumeric == 1);
[noinj]    = find(acl_pheno.InjuryNumeric == -1);
[APM]     = find(oa999_pheno.Disease == -1);
[OA]     = find(oa999_pheno.Disease == 1);


% FIGURE 6B
% Plot human PC scores in mouse space : Human age vs. PC4 colored by disease status
figure(1)
hold on
scatter( age(OA', 1), pc(OA', 4), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75);
scatter( age(APM', 1), pc(APM', 4),135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75);
xlabel('Age (years)'); 
ylabel('PC4'); 
title('ACLR OA vs. APM Human Scores in Mouse Space') ; 
legend('OA', 'APM'); 
xlim([25 85]); 
hold off

% FIGURE 3B
% Plot mouse PC1 vs. PC2 color coded by mouse injury status 
figure(2)
hold on
scatter( XS(inj', 1), XS(inj', 2), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75);
scatter( XS(noinj', 1), XS(noinj', 2),135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75);
xlabel('PC1'); 
ylabel('PC2'); 
title('ACLR Mouse PCA by Disease ');
legend('Injured', 'Uninjured'); 
axis equal;
hold off

