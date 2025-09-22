clear all; clc;

%load the model
load('model_ent.mat')   % COBRA version of Human-GEM v18 with entrez ID



% find index of exchange reactions (rxns) and related metabolites
rxns=printRxnFormula(model);
indexs_mtb=[]; %indexes of exchange reactions
mtb={}; %related metabolites
for i = 1:length(rxns)
    rxn=char(rxns(i));
    matches=strcmp('<=> ',rxn(end-3:end));
    if matches==1
    indexs_mtb=[indexs_mtb i];
    mtb{end+1}=rxn(1:end-6);
    end
end

% automatic way for finding index of exchange rxns
% [selExc, selUpt] = findExcRxns(model);
% indexs_mtb=find(selExc);

% close uptakes of all metabolites
model.lb(indexs_mtb)=0;


%find exchange rxns for metabolites of HAM's medium (reqired to maintain
%biomass production)
bio=importdata('rxn.xlsx'); %the file containing required rxns 
bio_rxns=bio(2:end,3); %required rxns 
%find indexes of required rxns
bio_index=[];
for i=1:length(bio_rxns)
    try
    bio_index(i)=find(strcmp(bio_rxns(i),model.rxns)); 
    catch
    end
end
%open and adjust boundary for uptakes of required metabolites
model.lb(bio_index)=-0.32/10;

% spesific constraints for some metabolite uptakes and release
model.lb(ismember(model.rxns,'MAR09034'))= -0.32;  % glugose uptake
model.lb(ismember(model.rxns,'MAR09048'))= -1.76;  % O2 uptake
model.lb(ismember(model.rxns,'MAR13082'))=0.0001;  % biomass release

% free uptake for HAM's medium metabolites that are not carbon source
model.lb(ismember(model.rxns,'MAR09076'))=-1000; %Fe HMR_9076
model.lb(ismember(model.rxns,'MAR09047'))=-1000;%H2O HMR_9047
model.lb(ismember(model.rxns,'MAR09072'))=-1000;%P HMR_9072
model.lb(ismember(model.rxns,'MAR09074'))=-1000;%sulfate HMR_9074

%include 3 extra metabolites that are known to be uptaken in brain
%metabolism
model.lb(ismember(model.rxns,'MAR09418'))=-0.32/10; % taurine
model.lb(ismember(model.rxns,'MAR09418'))=-0.32/10; % ornithnine
model.lb(ismember(model.rxns,'MAR09073'))=-0.32/10; % NH3 

%format the model for gurobi solver
model.A= [model.S]; % Stochiomatrix format for gurobi showen as A instead of S
model.A=sparse(model.A);
model.rhs=zeros(length(model.A(:,1)),1);
model.rxns=model.rxns; % Reaction names
model.obj=zeros(length(model.A(1,:)),1); % objective function
model.sense=['=']; 
model.vtype='C' ;

%save the contrained version of model to use in furher analysis
%save('model.constraint.mat','model')

%find the maximum release (also means max production) capacity of each
%external metabolite
max_mtb=[];
for i = 1:length(mtb)
    model2=model;
    model2.obj(indexs_mtb(i))=-1;
    mtb_rslt=gurobi(model2);
    max_mtb(i)=mtb_rslt.x(indexs_mtb(i));
end

% %save indexes and maximum capacities of external metabolites
% save('max_mtb.mat','max_mtb'); % save the maximum capacities
% save("indexs_mtb.mat","indexs_mtb"); %save indexes

%define metabolites that can be produced and released under given
%constraints and save thair names, indexes and max flux values

pro_ind=indexs_mtb(find(max_mtb>0));
pro_mtb=mtb(find(max_mtb>0));
max_pro=max_mtb(find(max_mtb>0));

save("pro_ind.mat","pro_ind");
save("pro_mtb.mat","pro_mtb");
save("max_pro.mat","max_pro");

% pro_mtb includes metabolites names with a coding system 
%translate metabolite names as human readable format and save
mtb={};
for i = 1:length(pro_mtb)
    mtb(i)=model.metNames(find(ismember(model.mets,pro_mtb(i))));
end
mtb=mtb'; 
save("met_names.mat","mtb"); %metabolites that can be produced and released