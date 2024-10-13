
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILES UTILIZED: TOAST_data.xlsx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Variable Workspace and Command Window
clc
clear

% Set up colors for heatmap
bwr = @(n)interp1([1 2 3], [0.512 0 0.128; 1 1 1; 1 1 1], linspace(1, 3, n), 'linear');
colormap(bwr(64));
map = [0 0 1
    1 0 0]; 

%% DMM oa007 (with control) heat maps

% FIGURE 3D
% GLM p-values for PC selection
DMM_ctrl_GLM = xlsread('TOAST_data.xlsx', 'DMM_w_control_GLMpvals'); 
figure(1)
yaxis = {'intercept', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13'} ; 
xaxis = {'mdl OA', 'mdl sex', 'mdl age'} ; 
h = heatmap(xaxis, yaxis, DMM_ctrl_GLM, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'DMM Control vs. OA GLM P-values by Model'; 

% FIGURE 5D
% TransComp-R p-values 
DMM_ctrl_TCR = xlsread('TOAST_data.xlsx', 'DMM_w_control_TCRpval'); 
figure(2)
yaxis = {'intercept', 'pc 3', 'pc 9', 'pc 11', 'pc 13', 'Age', 'Sex' ,'pc 3*Age', 'pc 9*Age', 'pc 11*Age', 'pc 13*Age', 'pc 3*Sex', 'pc 9*Sex', 'pc 11*Sex', 'pc 13*Sex'} ; 
xaxis = {'mdl pcs', 'mdl main', 'mdl age', 'mdl sex', 'mdl full'} ; 
h = heatmap(xaxis, yaxis, DMM_ctrl_TCR, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'DMM Control vs. OA Transcomp-R P-values by Model';

%% DMM oa999 (with APM) heat maps

% FIGURE 3F
% GLM p-values for PC selection
DMM_APM_GLM = xlsread('TOAST_data.xlsx', 'DMM_wAPM_GLMpvals'); 
figure(3)
yaxis = {'intercept', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13'} ; 
xaxis = {'mdl OA', 'mdl sex', 'mdl age'} ; 
h = heatmap(xaxis, yaxis, DMM_APM_GLM, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'DMM APM vs. OA GLM P-values by Model'; 

% FIGURE 5A
% TransComp-R p-values 
DMM_APM_TCR = xlsread('TOAST_data.xlsx', 'DMM_wAPM_TCRpval'); 
figure(4)
yaxis = {'intercept', 'pc 4', 'pc 8', 'pc 10', 'pc 13', 'Age', 'Sex' ,'pc 4*Age', 'pc 8*Age', 'pc 10*Age', 'pc 13*Age', 'pc 4*Sex', 'pc 8*Sex', 'pc 10*Sex', 'pc 13*Sex'} ; 
xaxis = {'mdl pcs', 'mdl main', 'mdl age', 'mdl sex', 'mdl full'} ; 
h = heatmap(xaxis, yaxis, DMM_APM_TCR, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'DMM APM vs. OA Transcomp-R P-values by Model';


%% ACLR oa007 (with control) heat maps

% FIGURE 3E
% GLM p-values for PC selection
ACLR_ctrl_GLM = xlsread('TOAST_data.xlsx', 'ACL_w_control_GLMpvals'); 
figure(5)
yaxis = {'intercept', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7'} ; 
xaxis = {'mdl OA', 'mdl sex', 'mdl age'} ; 
h = heatmap(xaxis, yaxis, ACLR_ctrl_GLM, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'ACLR Control vs. OA GLM P-values by Model'; 

% FIGURE 6C
% TransComp-R p-values 
ACLR_ctrl_TCR = xlsread('TOAST_data.xlsx', 'ACL_w_control_TCRpval'); 
figure(6)
yaxis = {'intercept', 'pc 1', 'pc 4', 'Age', 'Sex' ,'pc 1*Age', 'pc 4*Age', 'pc 1*Sex', 'pc 4*Sex'} ; 
xaxis = {'mdl pcs', 'mdl main', 'mdl age', 'mdl sex', 'mdl full'} ; 
h = heatmap(xaxis, yaxis, ACLR_ctrl_TCR, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'ACLR Control vs. OA Transcomp-R P-values by Model';


%% ACLR oa999 (with APM) heat maps

% FIGURE 3G
% GLM p-values for PC selection
ACLR_APM_GLM = xlsread('TOAST_data.xlsx', 'ACL_wAPM_GLMpvals'); 
figure(7)
yaxis = {'intercept', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7'} ; 
xaxis = {'mdl OA', 'mdl sex', 'mdl age'} ; 
h = heatmap(xaxis, yaxis, ACLR_APM_GLM, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'ACLR APM vs. OA GLM P-values by Model'; 

% FIGURE 6A
% TransComp-R p-values 
ACLR_APM_TCR = xlsread('TOAST_data.xlsx', 'ACL_wAPM_TCRpval'); 
figure(8)
yaxis = {'intercept', 'pc 3', 'pc 4', 'pc 6', 'Age', 'Sex' ,'pc 3*Age', 'pc 4*Age','pc 6*Age', 'pc 3*Sex', 'pc 4*Sex', 'pc 6*Sex'} ; 
xaxis = {'mdl pcs', 'mdl main', 'mdl age', 'mdl sex', 'mdl full'} ; 
h = heatmap(xaxis, yaxis, ACLR_APM_TCR, 'Colormap',colormap(bwr(100))) ; 
h.Title = 'ACLR APM vs. OA Transcomp-R P-values by Model';
