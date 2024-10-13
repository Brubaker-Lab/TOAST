
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

oa007_dataIn    = readtable('FinalData/GSE114007_human_log2_CPM.txt','Delimiter','\t','ReadRowNames',1);
oa007_pheno     = readtable('FinalData/GSE114007_human_phenotypes.txt','Delimiter','\t','ReadRowNames',0);

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
[~,h1,h2]           = intersect(oa007_dataIn.Properties.RowNames,oa007_dataIn.Properties.RowNames);
oa007_data          = oa007_dataIn(h1,:);

if ~isequal(oa007_data.Properties.RowNames,oa007_data.Properties.RowNames);
    disp('!!ERROR!! Human Genes Mismatch');end

[~,ih,im]           = intersect(oa007_data.Properties.RowNames,mouse_data.Properties.RowNames);
oa007_data          = oa007_data(ih,:);
mouse_data          = mouse_data(im,:);
if ~isequal(oa007_data.Properties.RowNames,mouse_data.Properties.RowNames);
    disp('!!ERROR!! ACL - OA007 Mismatch');end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACLR Mouse PCA and PC-R Model Building and Analysis
 
% Mouse PCA model for making the scores plots 
[~, ~, ~, ~, EXP]   = pca(zscore(mouse_data{:,:}')); 
numPC               = min(find(EXP < vThresh)) - 1; % Variance Threshold
[XL, XS, ~, ~, EXP] = pca(zscore(mouse_data{:,:}'),'NumComponents',numPC);

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
    exp_oa007(i,1)      =   XL(:,i)'*zscore(oa007_data{:,:}')'*zscore(oa007_data{:,:}')*XL(:,i)./...
        sum(diag(XL'*zscore(oa007_data{:,:}')'*zscore(oa007_data{:,:}')*XL));
end

% Variance Explained In Mouse and Human Model
vTable = splitvars(table([EXP(1:numPC), 100*exp_oa007]));
vTable.Properties.VariableNames = {'Mouse';'OA007'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merge Human Data Projections

human_data2      = oa007_data;
human_phen2      = oa007_pheno;

pc              = [zscore(human_data2{:,:}')*XL];
age             = [human_phen2.Age];
sex             = [human_phen2.Sex];
oa              = [human_phen2.Disease];
batch           = [ones(height(human_phen2),1)];

h_mdlTable                                          = splitvars(table(pc,age,sex,batch,oa));
h_mdlTable.Properties.VariableNames(numPC+1:end)    = {'Age';'Sex';'Batch';'Disease'};

% Human Models (GLM) used to pick PCs for the TransComp-R model
mdl_OA         = fitglm(h_mdlTable(:,[1:numPC end]));
mdl_set        = fitglm(h_mdlTable(:,[1:numPC end-1]));
mdl_sex        = fitglm(h_mdlTable(:,[1:numPC end-2]));
mdl_age        = fitglm(h_mdlTable(:,[1:numPC end-3]));

% Human TransComp-R Models (Combining Mouse PC's with Human Covariates) N = 58 human samples
pc_spec     = 'Disease ~ pc_1 + pc_4';
main_spec   = 'Disease ~ pc_1 + pc_4 + Age + Sex';
age_spec    = 'Disease ~ pc_1 + pc_4 + Age + Sex + pc_1*Age + pc_4*Age';
sex_spec    = 'Disease ~ pc_1 + pc_4 + Age + Sex + pc_1*Sex + pc_4*Sex';
full_spec   = 'Disease ~ pc_1 + pc_4 + Age + Sex + pc_1*Age + pc_4*Age + pc_1*Sex + pc_4*Sex';
simpleSpec = 'Disease ~ pc_1  + pc_4 + Age + Sex + pc_1*Sex + pc_1*Age + pc_4*Sex + pc_4*Age';

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

% Create lists of indices for various phenotypes
[inj]      = find(mouse_mdlTable.Disease == 1);
[noinj]    = find(mouse_mdlTable.Disease == -1);
[m_0day]     = find(mouse_mdlTable.Time == 0);
[m_1day]     = find(mouse_mdlTable.Time == 1);
[m_7day]     = find(mouse_mdlTable.Time == 7);
[m_14day]     = find(mouse_mdlTable.Time == 14);
[control]     = find(oa007_pheno.Disease == -1);
[OA]     = find(oa007_pheno.Disease == 1);

% FIGURE 3B
% Plot mouse PC1 vs. PC2 color coded by mouse injury status 
figure(1)
hold on
scatter(mouse_mdlTable.XS_1(noinj), mouse_mdlTable.XS_2(noinj), 135, 'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','b', 'LineWidth', 0.75) ; 
scatter(mouse_mdlTable.XS_1(inj), mouse_mdlTable.XS_2(inj), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor','r', 'LineWidth', 0.75); 
xlabel('XS1') ; 
ylabel('XS2') ; 
title('ACLR Mouse Linear Regression PCs by Injury ') ; 
legend('Uninjured', 'Injured');
xlim([-40 70])
hold off

% FIGURE 6D
% Plot human PC scores in mouse space : PC1 vs. Human Age colored by human disease status 
figure(2)
hold on
scatter( (age(control,1))', pc(control, 1),135, 'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor', 'b', 'LineWidth', 0.75);
scatter( (age(OA,1))', pc(OA, 1), 135,'MarkerEdgeColor',[ 0 0 0], 'MarkerFaceColor', 'r', 'LineWidth', 0.75);
xlabel('Age (years)') ; 
ylabel('PC1') ;  
title('ACLR OA vs. Control Human Scores in Mouse Space') ; 
xlim([10 80])
legend('Control', 'OA');
hold off

% create indices for PC1 boxplot groups 
num_M_noOA = 0; 
num_M_OA = 0; 
num_F_noOA = 0; 
num_F_OA = 0;

for i = 1: height(h_mdlTable)
    if h_mdlTable.Sex(i) == -1
        if h_mdlTable.Disease(i) == -1
            num_M_noOA = num_M_noOA + 1 ; 
            M_noOA(num_M_noOA, 1) = i;
        elseif h_mdlTable.Disease(i) == 1
            num_M_OA = num_M_OA + 1 ; 
            M_OA(num_M_OA, 1) = i;
        end
    elseif h_mdlTable.Sex(i) == 1
        if h_mdlTable.Disease(i) == -1
            num_F_noOA = num_F_noOA + 1 ; 
            F_noOA(num_F_noOA, 1) = i; 
        elseif h_mdlTable.Disease(i) == 1
            num_F_OA = num_F_OA + 1 ; 
            F_OA(num_F_OA, 1) = i;
        end
    end
end

% FIGURE 6E
% Boxplots of PC1 by phenotype
figure(3)
x = [h_mdlTable.pc_1(M_noOA',1); h_mdlTable.pc_1(M_OA',1); h_mdlTable.pc_1(F_noOA',1); h_mdlTable.pc_1(F_OA',1)];
g = [ones(size(h_mdlTable.pc_1(M_noOA',1))); 2*ones(size(h_mdlTable.pc_1(M_OA',1))); 3*ones(size(h_mdlTable.pc_1(F_noOA',1))); 4*ones(size(h_mdlTable.pc_1(F_OA',1)))];
boxplot(x, g, 'Labels',{'M_healthy','M_OA', 'F_healthy','F_OA'},'Whisker',1) ;
xlabel('Sex and Disease Status Groups') ;
ylabel('PC1') ; 
title('Boxplots of PC1 by Groups') ; 
ylim([-25 25])


%Mann Whitney test for PC1 boxplot
p1 = ranksum(h_mdlTable.pc_1(M_noOA',1), h_mdlTable.pc_1(M_OA',1)) ; 
p2 = ranksum(h_mdlTable.pc_1(M_noOA',1), h_mdlTable.pc_1(F_OA',1)) ; 
p3 = ranksum(h_mdlTable.pc_1(M_noOA',1), h_mdlTable.pc_1(F_noOA',1)) ; 
p4 = ranksum(h_mdlTable.pc_1(M_OA',1), h_mdlTable.pc_1(F_noOA',1)) ; 
p5 = ranksum(h_mdlTable.pc_1(M_OA',1), h_mdlTable.pc_1(F_OA',1)) ; 
p6 = ranksum(h_mdlTable.pc_1(F_OA',1), h_mdlTable.pc_1(F_noOA',1)) ; 

