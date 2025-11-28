clear
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_bacc.mat
load iAF_and_iYO_MCS.mat % your MCS file

run_all = 0;

%% Make the model
idx.bm = find(Com.c==1);
bio = find(Com.c==1);
idx.o2 = find(contains(Com.rxns,{'EX_o2_e'}));
idx.starch = find(contains(Com.rxns,'EX_starch_e'));
idx.glc  = find(ismember(Com.rxns,{'EX_glc__D_e_com'}));
Com.lb(idx.starch) = -1;
Com.lb(idx.glc) = 0;
% Setting the Cplex problem
%     Cplex_Com.Aeq=Com.S;
%     Cplex_Com.lb=Com.lb;
%     Cplex_Com.ub=Com.ub;
%     Cplex_Com.f=-Com.c;
%     Cplex_Com.beq=Com.b;
% % Run FBA with both organisms growth as objective
%     Cplex_Com.f(idx.bm(2)) = 0;
%     FBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
%     Cplex_Com.f(idx.bm(2)) = -1;
%     Cplex_Com.f(idx.bm(1)) = 0;
%     FBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
% % Set and run pFBA
%     Nogenes=cellfun(@isempty,Com.grRules);
%     Nogenes(find(ismember(Com.grRules,strcat(strcat('x(',num2str(find(contains(Com.genes,'shared')))),')')))) = 1; %Do not include the artificial genes
%     Genes=Nogenes==0;
%     Cplex_Com.lb(idx.bm) = FBA_org1(idx.bm);
%     Cplex_Com.f(:,1) = 0;
%     Cplex_Com.f(Genes,1)=1;
%     pFBA_org1=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
%     Cplex_Com.lb(idx.bm) = FBA_org2(idx.bm);
%     pFBA_org2=cplexlp(Cplex_Com.f,[],[],Cplex_Com.Aeq,Cplex_Com.beq,Cplex_Com.lb,Cplex_Com.ub);
%     
%     % Find all rxns that are always open in both cases
%     EX_rxns = find(startsWith(Com.rxns,'EX_'));
%     EX_Open_aero = [Com.rxns(EX_rxns(find(pFBA_org1(EX_rxns)~=0)));...
%         Com.rxns(EX_rxns(find(pFBA_org2(EX_rxns)~=0)))];
%     EX_Open_FBA = [Com.rxns(EX_rxns(find(FBA_org1(EX_rxns)~=0)));...
%         Com.rxns(EX_rxns(find(FBA_org2(EX_rxns)~=0)))];
%     Keep = [EX_Open_aero;...
%         Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org1'))));...
%         Com.rxns(find(ismember(Com.rxns,strrep(strrep(Com.rxns(find(startsWith(Com.rxns,'S_'))),'S_','EX_'),'_s','_e_org2'))))];
%     
%     % Close useless exchange
%     ub = max(Com.ub);
%     Com.lb(EX_rxns) = 0;
%     Com.ub(EX_rxns) = 0;
%     Com.lb(find(ismember(Com.rxns,Keep))) = Cplex_Com.lb(find(ismember(Com.rxns,Keep)));
%     Com.ub(find(ismember(Com.rxns,Keep))) = Cplex_Com.ub(find(ismember(Com.rxns,Keep)));
    Com_save = Com;
%% List all reactions in the MCS and perform analysis

  for i = 1:size(mcs,2)
    Solutions(1:size(find(ismember(mcs(:,i),-1)),1),i) = Com.rxns(find(ismember(mcs(:,i),-1)));
    mcs_size(i) = length(find(~ismember(mcs(:,i),0)));
  end
if run_all == 1;
  Good = 1;
  for i=1:size(mcs,2)
      Com = Com_save;
      Com.lb(mcs(:,i)==-1) = 0;
      Com.ub(mcs(:,i)==-1) = 0;
      Com.lb(find(endsWith(Com.rxns,'_s')))=0;
      Com.ub(find(endsWith(Com.rxns,'_s')))=0;
      FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
      Phenotype(1:2,i) = FBA(idx.bm);
      if Phenotype(1,i)>=0.05 && Phenotype(2,i) >=0.05
          Good(end+1) = i;
      end
  end
  Good(1) = [];
end 

%% finding minimal exchange phenotype
% Split shared exchange and import to find lowest flux strategy
MCS_chosen = 1;
Com = Com_save;
shared = find(contains(Com.rxns,'tr_s'));
ex_shared = find(endsWith(Com.rxns,'_s'));
ex_shared_name = Com.rxns(ex_shared);
ex_shared_S = Com.S(:,ex_shared);
Com.rxns(shared)=strcat(Com.rxns(shared),'_export');
Com.S(:,shared(end)+1:shared(end)+length(shared)) = -1*Com.S(:,shared);
Com.rxns(shared(end)+1:shared(end)+length(shared)) = strrep(Com.rxns(shared),'_export','_import');
Com.lb(shared(end)+1:shared(end)+length(shared)) = 0;
Com.ub(shared(end)+1:shared(end)+length(shared)) = 100;
Com.S(:,end+1:end+size(ex_shared_S,2)) = ex_shared_S;
Com.rxns(end+1:end+size(ex_shared_S,2)) = ex_shared_name;
Com.lb(end+1:end+size(ex_shared_S,2)) = 0;
Com.ub(end+1:end+size(ex_shared_S,2)) = 0;
Com.lb(find(contains(Com.rxns,'tr_s'))) = 0;
Com.c(end+1:length(Com.lb)) = 0;

% Add a sum of biomass variable
Com.mets(end+1) = {'sum_bio'};
Com.rxns(end+1) = {'bio_sum'};
Com.S(end+1,end+1) = 0;
Com.S(end,bio) = -1;
Com.S(end,end) = 1;
Com.b(end+1) = 0;
Com.lb(end+1) = 0;
Com.ub(end+1) = 100;
Com.c(end+1) = 0;
Com_save = Com; % this is the final matrix

% Calculate parsimonus maximal community growth rate
Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
Com.lb(find(endsWith(Com.rxns,'_s'))) = 0;
FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
Com.lb(find(ismember(Com.rxns,'bio_sum'))) = sum(FBA(bio));
Com.c(:) = 0;
Com.c(find(contains(Com.rxns,'tr_s'))) = -1;
Com.lb(find(contains(Com.rxns,'tr_s'))) = 0;
pFBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);

%% Phenotypic phase plane
remove = Com.rxns(mcs(:,MCS_chosen)==-1);
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
for i = 1:2
    Com.c(:) = 0;
    Com.c(bio(i)) = 1;
    FBA_self = cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    mu(:,i) = FBA_self(bio);
end
FBA_self(find(contains(Com.rxns,'EX_')));
step = 20;
bio1 = linspace(0,mu(1,1),step);
bio2 = linspace(0,mu(2,2),step);
growth_phpp = [];

for i = 1:step*2
    if i<=step
        Com = Com_save;
        Com.lb(find(ismember(Com.rxns,remove))) = 0;
        Com.ub(find(ismember(Com.rxns,remove))) = 0;
        Com.lb(bio(2)) = bio2(i);
        Com.ub(bio(2)) = bio2(i);
        Com.c(:) = 0;
        Com.c(bio(1)) = 1;
        FBAppp=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if isempty(FBAppp) == 0
            growth_phpp(i,1:2) = FBAppp(bio);
        else
            break
        end
    elseif i > step
        Com = Com_save;
        Com.lb(find(ismember(Com.rxns,remove))) = 0;
        Com.ub(find(ismember(Com.rxns,remove))) = 0;
        Com.lb(bio(1)) = bio1(i-step);
        Com.ub(bio(1)) = bio1(i-step);
        Com.c(:) = 0;
        Com.c(bio(2)) = 1;
        FBAppp=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if isempty(FBAppp) == 0
           growth_phpp(i,1:2) = FBAppp(bio);
        else
            break
        end
    end
end
