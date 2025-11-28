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
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_iAF1260_nar.mat

%% Make a list of all Ex_reactions that remains open
idx.bm = find(ismember(Com.c,1));
idx.o2 = find(contains(Com.rxns,{'EX_o2_e'}));
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
% Make a single product community export
if ~exist('cnan','var')
    startcna(1);
    wait(parfevalOnAll(@startcna,0,1));
end
ub = max(Com.ub);
Com.lb(EX_rxns) = 0;
Com.ub(EX_rxns) = 0;
Com.lb(find(ismember(Com.rxns,Keep))) = Cplex_Com.lb(find(ismember(Com.rxns,Keep)));
Com.ub(find(ismember(Com.rxns,Keep))) = Cplex_Com.ub(find(ismember(Com.rxns,Keep)));

% Make naringenins common
Com.S(find(contains(Com.mets,'narin_e_org1')),:) = ...
    Com.S(find(contains(Com.mets,'narin_e_org1')),:) + Com.S(find(contains(Com.mets,'narin_e_org2')),:);
Com.S(find(ismember(Com.mets,'narin_e_org2')),:) = 0;
Com.rxns(find(ismember(Com.rxns,'EX_narin_e_org1'))) = {'EX_narin_e_com'};
Com.mets(find(ismember(Com.mets,'narin_e_org1'))) = {'narin_e_com'};

Com.S(find(contains(Com.mets,'narin_e_org2')),:) = [];
Com.metNames(find(contains(Com.mets,'narin_e_org2'))) = [];
Com.metCharge(find(contains(Com.mets,'narin_e_org2'))) = [];
Com.metFormulas(find(contains(Com.mets,'narin_e_org2'))) = [];
Com.b(find(contains(Com.mets,'narin_e_org2'))) = [];
Com.mets(find(contains(Com.mets,'narin_e_org2'))) = [];

Com.S(:,find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.grRules(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.rxnNames(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.lb(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.ub(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.c(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.rev(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];
Com.rxns(find(contains(Com.rxns,'EX_narin_e_org2'))) = [];

% Close promiscuous
Com.ub(find(ismember(Com.rxns,'EX_narin_e_com'))) = ub;
Com.ub(find(ismember(Com.rxns,'CLcoa_org1'))) = 0;
Com.lb(find(ismember(Com.rxns,'CLcoa_org1'))) = 0;
Com.ub(find(ismember(Com.rxns,'TAL_org2'))) = 0;
Com.lb(find(ismember(Com.rxns,'TAL_org2'))) = 0;

% Make the CNA model
cnap = CNAcobra2cna(Com);

% Add gprRules from Cobra to CNA
%cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',Com.grRules);

cnap.reacMin(cnap.reacMin==-ub) = -inf;
cnap.reacMax(cnap.reacMax== ub) =  inf;
 
idx.glc  = find(ismember(cnap.reacID,{'EX_glc__D_e_com'}));
idx.nar  = find(ismember(cnap.reacID,{'EX_narin_e_com'}));
idx.atpm_org1 = find(ismember(cnap.reacID,{'ATPM_org1'}));
idx.atpm_org2 = find(ismember(cnap.reacID,{'ATPM_org2'}));
idx.bm = find(ismember(Com.c,1));
idx.o2 = find(contains(Com.rxns,{'EX_o2_e'}));

cnap.reacMax(idx.atpm_org1) = inf;
cnap.reacMax(idx.atpm_org2) = inf;

%% Uncomment to test feasibility 
 FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
 FBA(idx.bm)
 CNAplotPhasePlane(cnap,[idx.glc,idx.nar]);
%% Set up cMCS input

% Total glucose yield target (tested and working)
modules{1}.sense = 'target';
modules{1}.type = 'lin_constraints';
modules{1}.V = zeros(2,cnap.numr);
modules{1}.V(1,idx.nar) = 1;
modules{1}.V(2,idx.glc) = 1;
modules{1}.v = zeros(2,1);
modules{1}.v(1,1) = 0.1;
modules{1}.v(2,1) = -1;

modules{2}.sense = 'desired';
modules{2}.type = 'lin_constraints';
modules{2}.V = zeros(1,cnap.numr);
modules{2}.V(1,idx.bm(1)) = -1;
modules{2}.v = -0.1;

modules{3}.sense = 'desired';
modules{3}.type = 'lin_constraints';
modules{3}.V = zeros(1,cnap.numr);
modules{3}.V(1,idx.bm(2)) = -1;
modules{3}.v = -0.1;

% Promiscuous constraint
% modules{4}.sense = 'target';
% modules{4}.type = 'lin_constraints';
% modules{4}.V = zeros(2,cnap.numr);
% modules{4}.V(1,find(ismember(Com.rxns,'TAL_org1'))) = -1;
% modules{4}.V(2,find(ismember(Com.rxns,'CLcoa_org1'))) = -1;
% modules{4}.v = zeros(2,1);
% modules{4}.v(1,1) = -0.0001;
% modules{4}.v(2,1) = -0.0001;
% 
% Promiscuous constraint
% modules{5}.sense = 'target';
% modules{5}.type = 'lin_constraints';
% modules{5}.V = zeros(2,cnap.numr);
% modules{5}.V(1,find(ismember(Com.rxns,'TAL_org2'))) = -1;
% modules{5}.V(2,find(ismember(Com.rxns,'CLcoa_org2'))) = -1;
% modules{5}.v = zeros(2,1);
% modules{5}.v(1,1) = -0.0001;
% modules{5}.v(2,1) = -0.0001;

%% Rxns that cannot be ko'ed
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
  
  % Set all lethal to nan
     rkoCost = ones(cnap.numr,1);
     Genes = ~cellfun(@isempty,Com.grRules);
     rkoCost(find(ismember(Genes,0))) = nan; % nan means its excluded from the possible solution
     rkoCost(find(ismember(Com.grRules,'s0001'))) = nan;
     rkoCost(find(contains(Com.rxns,'tr_s_'))) = nan;
     rkoCost(find(contains(Com.rxns,'CO2tex_'))) = nan;
     rkoCost(find(contains(Com.rxns,'Htex_'))) = nan;
     rkoCost(lethal_rxn) = nan;
    
  % Make the nar reactions feasible
    rkiCost = nan(cnap.numr,1);
%     nar_rxns = {'TAL_org1','TAL_org2',...
%         'Pa4CL_org1','Pa4CL_org2','CHS_org1','CHS_org2',...
%         'CHI_org1','CHI_org2'};
%     cnap.reacMin(find(ismember(Com.rxns,nar_rxns))) = 0;
%     cnap.reacMax(find(ismember(Com.rxns,nar_rxns))) = 0;
%     rkiCost(find(ismember(Com.rxns,nar_rxns))) = 1;
%     rkiCost(find(ismember(Com.rxns,trs_nar))) = 1;


% rkiCost(find(contains(Com.rxns,'_s_org'))) = 1;
% specifying gene knockout costs
% Generate list of genes as a template to define k/o genes
% [~,~,genes,gpr_rules] = CNAgenerateGPRrules(cnap);
% gkoCost = ones(numel(genes),1);
% gkiCost = nan(numel(genes),1);
% gkiCost(find(ismember(genes,strcat('x',num2str(find(contains(Com.genes,'shared'))))))) = 0;

maxSolutions = 200;
maxCost = 25;
options.milp_solver = 'cplex';
options.milp_time_limit  = 43200;
%options.milp_bigM = 100000; %Try true, options.seed (change could help)
options.mcs_search_mode = 1;
verbose = 4;

%% Run cMCS and save outcome 

[mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
    CNAMCSEnumerator3(cnap,modules,rkoCost,rkiCost,maxSolutions,maxCost,options,verbose);

if isempty(mcs) == 0
    mcs(isnan(mcs)) = 0;
    cnap.reacMin(mcs(:,end)==-1) = 0;
    cnap.reacMax(mcs(:,end)==-1) = 0;
    CNAplotPhasePlane(cnap,[idx.bm;idx.nar]);
    FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
    FBA(idx.bm)
    FBA(idx.nar)
    cd /nfs/homes/tremb133/CommunityStuff/ComMCS
    save ('Jan_NAR_MCS_core_real.mat','mcs')
else
    'no solution'
end