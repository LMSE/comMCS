clear;
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
addpath('/nfs/homes/tremb133/CommunityStuff/ComMCS');
lethal = readtable('Lethal_genes.csv');
lethal = table2array(lethal);
% Cplex Path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux')
% CNA path
addpath('/nfs/homes/tremb133/CellNetAnalyzerNEWER')
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_iAF.mat

%% Make a list of all Ex_reactions to be close open
idx.bm = find(Com.c==1);
idx.o2 = find(contains(Com.rxns,{'EX_o2_e'}));
idx.glc  = find(ismember(Com.rxns,{'EX_glc__D_e_com'}));
% Setting the Cplex problem
    Cplex_Com.Aeq=Com.S;
    Cplex_Com.lb=Com.lb;
    Cplex_Com.ub=Com.ub;
    Cplex_Com.f=-Com.c;
    Cplex_Com.beq=Com.b;
% Run FBA with both organisms growth as objective
    Cplex_Com.f(idx.bm(2)) = 0;
    FBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    Cplex_Com.f(idx.bm(2)) = -1;
    Cplex_Com.f(idx.bm(1)) = 0;
    FBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
% Set and run pFBA
    Nogenes=cellfun(@isempty,Com.grRules);
    Nogenes(find(ismember(Com.grRules,strcat(strcat('x(',num2str(find(contains(Com.genes,'shared')))),')')))) = 1; %Do not include the artificial genes
    Genes=Nogenes==0;
    Cplex_Com.lb(idx.bm) = FBA_org1(idx.bm);
    Cplex_Com.f(:,1) = 0;
    Cplex_Com.f(Genes,1)=1;
    pFBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    Cplex_Com.lb(idx.bm) = FBA_org2(idx.bm);
    pFBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
    
    % Find all rxns that are always open in both cases
    EX_rxns = find(startsWith(Com.rxns,'EX_'));
    EX_Open_aero = [Com.rxns(EX_rxns(find(pFBA_org1(EX_rxns)~=0)));...
        Com.rxns(EX_rxns(find(pFBA_org2(EX_rxns)~=0)))];
    EX_Open_FBA = [Com.rxns(EX_rxns(find(FBA_org1(EX_rxns)~=0)));...
        Com.rxns(EX_rxns(find(FBA_org2(EX_rxns)~=0)))];
    Keep = [EX_Open_aero;...
        Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org1'))));...
        Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org2'))))];

%% Make the model

if ~exist('cnan','var')
    startcna(1);
    wait(parfevalOnAll(@startcna,0,1));
end
ub = max(Com.ub);
Com.lb(EX_rxns) = 0;
Com.ub(EX_rxns) = 0;
Com.lb(find(ismember(Com.rxns,Keep))) = Cplex_Com.lb(find(ismember(Com.rxns,Keep)));
Com.ub(find(ismember(Com.rxns,Keep))) = Cplex_Com.ub(find(ismember(Com.rxns,Keep)));

% Make the CNA model
cnap = CNAcobra2cna(Com);

% Add gprRules from Cobra to CNA
cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',Com.grRules);

cnap.reacMin(cnap.reacMin==-ub) = -inf;
cnap.reacMax(cnap.reacMax== ub) =  inf;

idx.atpm_org1 = find(ismember(cnap.reacID,{'ATPM_org1'}));
idx.atpm_org2 = find(ismember(cnap.reacID,{'ATPM_org2'}));

cnap.reacMax(idx.atpm_org1) = inf;
cnap.reacMax(idx.atpm_org2) = inf;

%% Uncomment to test feasibility 
%  FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
%  FBA(idx.bm)
%  CNAplotPhasePlane(cnap,idx.bm)

%% Set up cMCS input
modules{1}.sense = 'target';
modules{1}.type = 'lin_constraints';
modules{1}.V = zeros(2,cnap.numr);
modules{1}.V(1,idx.bm(1)) = 1;
modules{1}.V(2,idx.bm(2)) = -1;
modules{1}.v = zeros(2,1);
modules{1}.v(2,1) = -0.01;

modules{2}.sense = 'target';
modules{2}.type = 'lin_constraints';
modules{2}.V = zeros(2,cnap.numr);
modules{2}.V(1,idx.bm(2)) = 1;
modules{2}.V(2,idx.bm(1)) = -1;
modules{2}.v = zeros(2,1);
modules{2}.v(2,1) = -0.01;

modules{3}.sense = 'desired';
modules{3}.type = 'lin_constraints';
modules{3}.V = zeros(1,cnap.numr);
modules{3}.V(1,idx.bm(1)) = -1;
modules{3}.v = zeros(1,1);
modules{3}.v(:,:) = -0.0125;

modules{4}.sense = 'desired';
modules{4}.type = 'lin_constraints';
modules{4}.V = zeros(1,cnap.numr);
modules{4}.V(1,idx.bm(2)) = -1;
modules{4}.v = zeros(1,1);
modules{4}.v(:,:) = -0.0125;
%% Lethal Rxns should not be feasible ko
% remove experimentally found lethal kos
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
  
  % Set atll lethal to nan
    rkoCost = ones(cnap.numr,1);
    Genes = ~cellfun(@isempty,Com.grRules);
    rkoCost(find(ismember(Genes,0))) = nan; % nan means its excluded from the possible solution
    rkoCost(find(contains(Com.grRules,'s0001'))) = nan;
    rkoCost(find(contains(Com.rxns,'tr_s_'))) = nan;
    rkoCost(find(contains(Com.rxns,'CO2tex_'))) = nan;
    rkoCost(find(contains(Com.rxns,'Htex_'))) = nan;
    rkoCost(lethal_rxn) = nan;
    
    rkiCost = nan(cnap.numr,1); % Here no ki
% rkiCost(find(contains(Com.rxns,'_s_org'))) = 1;
% specifying gene knockout costs
% Generate list of genes as a template to define k/o genes
% [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
% gkoCost = ones(numel(genes),1);
% gkiCost = nan(numel(genes),1);
% gkiCost(find(ismember(genes,strcat('x',num2str(find(contains(Com.genes,'shared'))))))) = 0;

maxSolutions = 1000000;
maxCost = 50;
options.milp_solver = 'cplex';
options.milp_time_limit  = 3600*24;
options.milp_bigM = 1000000; %Try true, options.seed (change could help)
options.mcs_search_mode = 1;
verbose = 1;

%% Run cMCS and save outcome 

[mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
    CNAMCSEnumerator3(cnap,modules,rkoCost,rkiCost,maxSolutions,maxCost,options,verbose);

% [mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
%     CNAgeneMCSEnumerator3(cnap,modules,{},{},maxSolutions,maxCost,gkoCost,gkiCost,{},options,verbose);
Com.keep = Keep;
cd /nfs/homes/tremb133/CommunityStuff/ComMCS
save ('October20_iAF_MCS.mat','mcs');
%save ('MCS_Com_iAF.mat','Com')  
mcs_test = mcs(:,1);
cnap.reacMin(mcs_test == -1) = 0;
cnap.reacMax(mcs_test == -1) = 0;
CNAplotPhasePlane(cnap,idx.bm)
cnap.objFunc(idx.bm(1)) = 0;
FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
FBA(idx.bm)