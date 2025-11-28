clear;
cd '/nfs/homes/tremb133/CommunityStuff/Ratio/model';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load MCS_Com_iAF
load October_iAF_MCS
load eco_iAF1260 % individual models
load Gene_mcs_iAF_keio_Oct.mat
load 'MCS_valid_iAF_keio_Oct.mat'
lethal = readtable('Lethal_genes.csv');
lethal = table2array(lethal);
model = iAF1260;
model_save = model;
bio = find(ismember(model.c,1));
Analysis.Lethal_org1 = [];
Analysis.Lethal_org2 = [];
Analysis.Growth_after_MCS_org1 = [];
Analysis.Growth_after_MCS_org2 = [];
Analysis.Strategy_org1 = [];
Analysis.Strategy_org2 = [];

%% Find what each individual strain is dependant of
MetEx = {'ade','ala__L','ala__D','akg','arg__L','arg__D','asn__L',...
     'asn__D','asp__L','asp__D','btn','fol','glu__L','glu__D',...
     'gln__L','gln__D','gly','leu__L','leu__D','ile__L','met__L',...
     'met__D','pnto__R','phe__L','phe__D','purine','pyrimid',...
     'ribflv','ser__L','ser__D','thr__L','tyr__L','thym','trp__L','ura'};

%real_removal = Solutions;
Ex_loc = find(ismember(model.rxns, strcat('EX_',[strcat(MetEx,'_e')])));
Analysis.Mets = model.rxns(Ex_loc)';

% remove lethal rxns
lethal_rxn = 0;
  for i = 1:length(lethal)
      rxn_w_gene = find(contains(Com.grRules,lethal(i)));
      for j = 1:length(rxn_w_gene)
          if contains(Com.grRules(rxn_w_gene(j)),'and') == 1
              continue
          else
              lethal_rxn(end+1) = rxn_w_gene(j);
          end
      end     
  end
lethal_rxn(1) = [];
lethal_rxn = unique(lethal_rxn);

for i = 1:size(real_removal,2)
    if any(ismember(Com.rxns(lethal_rxn),real_removal(find(~cellfun(@isempty,real_removal(:,i))),i)))
        works(i) = 0;
    else
        works(i) = 1;
    end
end
valid = find(works==1);    
% =========== Organism 1 =========== %
for i=1:size(real_removal,2)
   if any(valid==i)
       remove = real_removal(:,i);
       rxns = remove(~cellfun(@isempty,real_removal(:,i)));
       org1_rmv = rxns(find(contains(rxns,'_org1')));
       org1_rmv = strrep(org1_rmv,'_org1','');
       model = model_save;
       model.lb(find(ismember(model.rxns,org1_rmv)))=0;
       model.ub(find(ismember(model.rxns,org1_rmv)))=0;
       FBA=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);
       Analysis.Growth_after_MCS_org1(end+1) = FBA(bio);
       Analysis.Strategy_org1(end+1) = i;
       if FBA(bio) > 0
          Analysis.Lethal_org1(end+1) = 0;
       else
          Analysis.Lethal_org1(end+1) = 1;
       end
   end
end

% =========== Organism 2 =========== %
for i=1:size(real_removal,2)
   if any(valid==i)
        remove = real_removal(:,i);
        rxns = remove(~cellfun(@isempty,real_removal(:,i)));
        org2_rmv = rxns(find(contains(rxns,'_org2')));
        org2_rmv = strrep(org2_rmv,'_org2','');
        model = model_save;
        model.lb(find(ismember(model.rxns,org2_rmv)))=0;
        model.ub(find(ismember(model.rxns,org2_rmv)))=0;
        FBA=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);
        Analysis.Growth_after_MCS_org2(end+1) = FBA(bio);
        Analysis.Strategy_org2(end+1) = i;
        if FBA(bio) > 0
            Analysis.Lethal_org2(end+1) = 0;
        else
            Analysis.Lethal_org2(end+1) = 1;
        end
   end
end
clear i
to_test = [];
for i = 1:length(Analysis.Lethal_org1)
    if Analysis.Lethal_org1(i) == 1 && Analysis.Lethal_org2(i) == 1
        to_test(end+1) = Analysis.Strategy_org1(i);
    end
end
cd '/nfs/homes/tremb133/CommunityStuff/ComMCS';
save ('To_test_iAF_keio_Oct.mat','to_test');