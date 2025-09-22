clear all; clc;
tic

%default FCs
default_fc=ones(12995,1);
load('indexs_mtb.mat'); %indexes of exchange reactions
default_fc(indexs_mtb)=2;

%load required data
load('max_pro.mat') %max flux for metabolites that are produced and released
load('pro_ind.mat') %indexes of metabolites that are produced and released
load('pro_mtb.mat') %metabolites that are produced and released
load('met_names.mat') %readable format of metabolite names
load('Level_Data_GSE20292.mat') %example data named as Level_Data 
%it is created in map_rxn_level.m (you can use your own data)
%Level_Data includes:
%Level_Data.FC that include transcriptome based FC values for each sample
%in the responsible dataset

% Level_Data.FC(find(Level_Data.FC<0.1))=0.1;
% Level_Data.FC(find(Level_Data.FC>10))=10;

for i=1:size(Level_Data.FC,2) %for each sample of dataset

    fc=Level_Data.FC(:,i);

    for j=1:2
    
        load('model_constraint.mat') %load model constrained in the previous part
        % model.exch=ones(12995,1);
        % model.exch(indexs_mtb)=0;
        %modele yeni bir değişken ekleyip exchage rxnlar için 0 diğerleri için
        %1 girirerek daha sonra sonuçlarda exchangleri toplamdan exclude
        %edebiliz
        
        %enter fc sets, * or / depending on condition
        if j==1
            model.fc=default_fc.*fc; %for control case
        else
            model.fc=default_fc./fc; %for disease case
        end

        %irreversible form of model
        [modelIrrev,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);
        model=modelIrrev;
        
        %format for gurobi
       model.A= [model.S]; % Stochiomatrix format for gurobi showen as A instead of S
       model.A=sparse(model.A);
       model.rhs=zeros(length(model.A(:,1)),1);
       model.rxns=model.rxns; % Reaction names
       model.obj=zeros(length(model.A(1,:)),1); % objective function
       model.sense=['=']; 
       model.vtype='C' ;

       
       % org_scr=[];
       % int_scr=[];
       % rxn_act4=[];
       rxn_act5=[];
       %parpool(16); %delete(gcp('nocreate'));
       %TAMBOOR: minimization of weighted reactions
       %FCs are used as reactions wieghts in objective function
       %make it for each metabolites that can be produced and released 
       parfor k = 1:length(pro_ind)
            model2=model;
            model2.lb(pro_ind(k))=max_pro(k)*0.9;
            model2.obj=model2.fc;
            rslt=gurobi(model2);             
            % int_scr(k)=sum(rslt.x(find(model2.exch==1))); 
            % rxn_act4(k)=length(find(rslt.x>0.0001));
            rxn_act5(k)=length(find(rslt.x>0.00001)); %number of active rxn
            % orj_scr(k)=rslt.objval;
        end
        % orj(:,j)=orj_scr';
        % ints(:,j)=int_scr';
        % act4(:,j)=rxn_act4';
        act5(:,j)=rxn_act5';
    end

    ScoreData(i).DataID=Level_Data.DataID;
    ScoreData(i).metNames=mtb;

    % %for orjinal
    % scores=orj(1:end,:);
    % raw_score=(scores(:,1) - scores(:,2))./(scores(:,1) + scores(:,2));
    % ort=mean(raw_score);
    % st_d=std(raw_score);
    % z_score=(raw_score - ort)./st_d;
    % ScoreData(i).Org_timbrScore(:,1)=scores(:,1);
    % ScoreData(i).Org_timbrScore(:,2)=scores(:,2);
    % ScoreData(i).Org_timbrScore(:,3)=z_score;
    % 
    % %for internal sum
    % scores=ints(1:end,:);
    % raw_score=(scores(:,1) - scores(:,2))./(scores(:,1) + scores(:,2));
    % ort=mean(raw_score);
    % st_d=std(raw_score);
    % z_score=(raw_score - ort)./st_d;
    % ScoreData(i).SumInt_timbrScore(:,1)=scores(:,1);
    % ScoreData(i).SumInt_timbrScore(:,2)=scores(:,2);
    % ScoreData(i).SumInt_timbrScore(:,3)=z_score;
    % 
    % %for act4
    % scores=act4(1:end,:);
    % raw_score=(scores(:,1) - scores(:,2))./(scores(:,1) + scores(:,2));
    % ort=mean(raw_score);
    % st_d=std(raw_score);
    % z_score=(raw_score - ort)./st_d;
    % ScoreData(i).act4_timbrScore(:,1)=scores(:,1);
    % ScoreData(i).act4_timbrScore(:,2)=scores(:,2);
    % ScoreData(i).act4_timbrScore(:,3)=z_score;

    %calculate production scores and turn them z-scores
    scores=act5(1:end,:);
    raw_score=(scores(:,1) - scores(:,2))./(scores(:,1) + scores(:,2));
    ort=mean(raw_score);
    st_d=std(raw_score);
    z_score=(raw_score - ort)./st_d;
    ScoreData(i).act5_timbrScore(:,1)=scores(:,1);
    ScoreData(i).act5_timbrScore(:,2)=scores(:,2);
    ScoreData(i).act5_timbrScore(:,3)=z_score;
end
toc

save("ScoreData_GSE20292.mat","ScoreData") %save the results for that dataset