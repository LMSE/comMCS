clear
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_iMLcore_nar.mat
load Jan_NAR_MCS_core_real.mat % your MCS file

run_all = 0;

%% Make the model
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
    
    Com.ub(find(ismember(Com.rxns,'CLcoa_org1'))) = 0;
    Com.lb(find(ismember(Com.rxns,'CLcoa_org1'))) = 0;
    Com.ub(find(ismember(Com.rxns,'TAL_org2'))) = 0;
    Com.lb(find(ismember(Com.rxns,'TAL_org2'))) = 0;

    Com.ub(find(ismember(Com.rxns,'EX_narin_e_com'))) = 100;
    % Find components location
    bio = find(ismember(Com.c,1));
    nar = find(ismember(Com.rxns,'EX_narin_e_com'));
    Com_save = Com;
    
%% List all reactions in the MCS and perform analysis

  for i = 1:size(mcs,2)
    Solutions(1:size(find(ismember(mcs(:,i),-1)),1),i) = Com.rxns(find(ismember(mcs(:,i),-1)));
    mcs_size(i) = length(find(~ismember(mcs(:,i),0)));
  end
nar_ex = {'EX_coumaci__p_s','EX_narin_s','EX_tyr__L_s'};
nar_rxns = {'TAL_','CLcoa_','CHS_','CHI_'};
if run_all == 1
  Good = 1;
  %for i=1:8
  for i=1:size(mcs,2)
      Com = Com_save;
      Com.lb(mcs(:,i)==-1) = 0;
      Com.ub(mcs(:,i)==-1) = 0;
      Com.lb(find(endsWith(Com.rxns,'_s')))=0;
      Com.ub(find(endsWith(Com.rxns,'_s')))=0;
      FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
      Phenotype(1:2,i) = FBA(bio);
      Phenotype(3,i) = FBA(nar);
      if Phenotype(1,i)>=0.001 && Phenotype(2,i) >=0.001 && Phenotype(3,i) > 0
         Good(end+1) = i;
         for j=1:length(nar_rxns)
            if sum(FBA(find(contains(Com.rxns,nar_rxns(j))))>0)<3
              org(j,1) = find(FBA(find(contains(Com.rxns,nar_rxns(j))))>0);
            else
              org(j,1) = 3;
            end
         end
      else
          org(1:length(nar_rxns),1)=0;
      end
      nar_org(:,i) = org;
  end
  Good(1) = [];
end 
if run_all ==1
    for i = 1:length(Good)
        num_strat(i) = length(unique(strrep(strrep(Solutions(~cellfun(@isempty,Solutions(:,i)),i),'_org1',''),'_org2','')));
    end
end
%% Phenotypic phase plane
% Product
MCS_chosen = 779;
glu = find(ismember(Com.rxns,'EX_glc__D_e_com'));
Com = Com_save;
Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;

for i = 1:2
    Com.c(:) = 0;
    Com.c(bio(i)) = 1;
    FBA_self = cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    mu(:,i) = FBA_self(bio);
end
bio1 = linspace(0,mu(1,1),10);
bio2 = linspace(0,mu(2,2),10);
glu_up = linspace(-10,0,10);
Com.c(:) = 0;
Com.c(nar) = 1;
FBA_nar = cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
nar_up = linspace(0,FBA_nar(nar),10);
nar_production = [];
gluc_uptake = [];
strain1 = [];
strain2 = [];
for i = 1:20
    if i<=10
        Com = Com_save;
        Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
        Com.lb(find(endsWith(Com.rxns,'_s'))) = 0;
        Com.lb(glu) = glu_up(i);
        Com.c(:) = 0;
        Com.c(nar) = 1;
        FBAppp=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if isempty(FBAppp) == 0
            nar_production(i,1) = FBAppp(nar);
            gluc_uptake(i,1) = FBAppp(glu);
        else
            break
        end
    elseif i > 10
        Com = Com_save;
        Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
        Com.lb(find(endsWith(Com.rxns,'_s'))) = 0;
        Com.lb(nar) = nar_up(i-10);
        Com.c(:) = 0;
        Com.c(glu) = -1;
        FBAppp=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if isempty(FBAppp) == 0
            nar_production(i,1) = FBAppp(nar);
            gluc_uptake(i,1) = FBAppp(glu);
        else
            break
        end
    end
end

% biomass
Com = Com_save;
Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
for i = 1:2
    Com.c(:) = 0;
    Com.c(bio(i)) = 1;
    FBA_self = cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    mu(:,i) = FBA_self(bio);
end
FBA_self(find(contains(Com.rxns,'EX_')));
step = 10;
bio1 = linspace(0,mu(1,1),step);
bio2 = linspace(0,mu(2,2),step);
growth_phpp = [];

for i = 1:step*2
    if i<=step
        Com = Com_save;
        Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
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
        Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
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
%% finding minimal exchange phenotype
% Split shared exchange and import to find lowest flux strategy
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

%% dFBA
met_list = upper(erase(erase(Com.rxns(find(endsWith(Com.rxns,'_s'))),'EX_'),'_s'));
% optimal scenario
    % Here initially, for the first hour we will assume they have enough
    % things to import from the environment, then we calculate mu, allowing
    % us to calculate the X, which will be used to calculate the ub on imporation.
    % set-up initial parameters
    V=1;
    Yxs = 0.85;
    d = 0.1/3600;
    X0 = [0.1;0.1]; % initial biomass (g-DW/L)
    S0 = 10/180*1000; % inital glucose concentration (mmol/L)
    dt = 0.1; % space between iteration
    t0= 0;
    kd = 0.01 ;% decay
    km = 0.0001; % m-m kinetics
    Com = Com_save;
    Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
    Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
    Com.ub(find(ismember(Com.rxns,nar_ex))) = 0;
    % Reset variables and initialize
    X=0;
    P=0;
    S=0;
    i=1;
    S(i) = S0;
    X(1:2,i) = X0;
    met_conc(1:length(met_list),1) = 0;
    Glu_Vmax = Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))); % Depends on model
    % Begin dFBA
        while isempty(FBA)==0
            Glu_lb = round(-Glu_Vmax*(S(i)/(S0+4/180/1000)),3); % 4 mg/L Monod 1946
            Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))) = -1*Glu_lb;
            FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
            if sum(FBA(bio)) <= 0.0001 | i > 1000
                break
            end
            X(1:2,i+1) = X(:,i).*exp((FBA(bio)-kd)*dt);
            S(i+1) = S(i)-Glu_lb*sum(X(:,i))*dt;
            P(i+1) = P(i)+FBA(find(ismember(Com.rxns,'EX_narin_e_com')))*sum(X(:,i))*dt;
            % Set-up organism ratio constraint, except for shared part,
            % growth and ATPM
            Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.ub(bio) = 100;
            Com.lb(bio) = 0;
            Com.ub(find(contains(Com.rxns,'ATPM'))) = 100;
            % shared is the minimal scenario captured by pFBA
%             Com.ub(find(contains(Com.rxns,'tr_s'))) = abs(pFBA(find(contains(Com.rxns,'tr_s'))));

            % Neglect the first 10 hour, as we assume media is enough
            if i <= 100
                importation(1:length(met_list),1) = 0;
                importation(1:length(met_list),2) = 0;   
                importation(9,1) = abs(pFBA(find(contains(Com.rxns,'COUMACI__Ptr_s_org1_import'))));
                importation(9,2) = abs(pFBA(find(contains(Com.rxns,'COUMACI__Ptr_s_org2_import'))));
            else
                importation(:,1) = abs(pFBA(find(contains(Com.rxns,'tr_s_org1_import'))));
                importation(:,2) = abs(pFBA(find(contains(Com.rxns,'tr_s_org2_import'))));
            end   
            
            production(:,1) = abs(pFBA(find(contains(Com.rxns,'tr_s_org1_export'))));
            production(:,2) = abs(pFBA(find(contains(Com.rxns,'tr_s_org2_export'))));
            % Perform one iteration 
                for j = 1:length(met_list)
                    dx(j,i) = production(j,1)*X(1,i) + production(j,2)*X(2,i)-...
                        importation(j,1)*X(1,i) - importation(j,2)*X(2,i);
                    conc(j) = met_conc(j,i) + dx(j,i)*dt;
                    if conc(j) < 0
                       conc(j) = 0;
                    end
                    if i>=100
                        loc = find(contains(Com.rxns,strcat(met_list(j),'tr_s'))); % location of met
%                         Com.ub(loc(3)) = abs(pFBA(loc(3))*(conc(j)/(conc(j)+km))*X(1,i+1)/sum(X(:,i+1)));
%                         Com.ub(loc(4)) = abs(pFBA(loc(4))*(conc(j)/(conc(j)+km))*X(2,i+1)/sum(X(:,i+1)));
                    elseif i>=10
                        if j == 9
                            loc = find(contains(Com.rxns,strcat(met_list(j),'tr_s')));
%                             Com.ub(loc(3)) = abs(pFBA(loc(3))*(conc(j)/(conc(j)+km))*X(1,i+1)/sum(X(:,i+1)));
%                             Com.ub(loc(4)) = abs(pFBA(loc(4))*(conc(j)/(conc(j)+km))*X(2,i+1)/sum(X(:,i+1)));
                        end
                    end
                end
            Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
            Com.lb(ismember(mcs(:,MCS_chosen),-1))=0;
            Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
            met_conc(:,i+1) = conc;
            i=i+1;
        end
gluc_gl = S'*0.180;
nar_gl = P'*0.272;
X_graph = X';