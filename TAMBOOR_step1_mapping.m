clear all; clc;

load('model_ent.mat'); % COBRA version of Human-GEM v18 with entrez ID


load('SNData_All2.mat'); %(SNDataAll) includes normalized expression data of each samples
%you can use your own data
SNData=SNData_All;
%SNData includes:
%SNData.DataID -> ID for dataset
%SNData.GeneID -> Entrez Gene IDs for related dataset
%SNData.Control -> Normalized gene expression values for control samples
%SNData.PD -> Normalized gene expression values for PD samples

%for i=1:8 % you can apply algortihm for all data with loop if you want

%make it separately for each dataset
%for example data for GSE20292 (it is in 4th order in mat file (SNData))
%for control case map mean expression value
exp_data=[];
exp_data.gene=SNData(4).GeneID;
exp_data.value= table2array(mean(SNData(4).Control,2)); 
[exp_rxns parsedGPR] = mapExpressionToReactions(model, exp_data);  %Expression data will be mapped to reactions of model.

% Find NaN values using isnan function
nanIndices = isnan(exp_rxns);
% Replace NaN values with 1 using logical indexing
exp_rxns(nanIndices) = 1;
 %mapped data (reaction scores) are saved in new variable
Level_Data.DataID = SNData(4).DataID; 
Level_Data.C = exp_rxns;

%for PD case
exp_data=[];
exp_data.gene=SNData(4).GeneID;
for i =1:size(SNData(4).PD,2) %for each samples of PD case
    exp_data.value= table2array(SNData(4).PD,2); 
    exp_data.value= exp_data.value(:,i);
    [exp_rxns parsedGPR] = mapExpressionToReactions(model, exp_data);  %Expression data willbe mapped to reactions of model.
   
    % Find NaN values using isnan function
    nanIndices = isnan(exp_rxns);
    % Replace NaN values with 1 using logical indexing
    exp_rxns(nanIndices) = 1;
    Level_Data.DataID = SNData(4).DataID;
    Level_Data.PD(:,i) = exp_rxns;
    %calculate FC for each PD samples
    Level_Data.FC(:,i) = Level_Data.PD(:,i)./Level_Data.C;
end

%save the reaction scores and FC 
save("Level_Data_GSE20292.mat","Level_Data");








