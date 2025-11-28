clear;
cd '/nfs/homes/tremb133/CommunityStuff/Ratio/model';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load MCS_Com_iMLcore_OLD.mat
load ComStrategy_4_Core
bio = find(Com.c==1);

%% List all reactions in MCS

for i = 1:size(mcs,2)
    Solutions(1:size(find(ismember(mcs(:,i),-1)),1),i) = Com.rxns(find(ismember(mcs(:,i),-1)));
    mcs_size(i) = length(find(~ismember(mcs(:,i),0)));
end
MCS_chosen =49;
remove = Solutions(~cellfun(@isempty, Solutions(:,MCS_chosen)),MCS_chosen);
%% Calculate minimal parsimonious exchanges
% Split shared exchange and import
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
exchange_close = {'EX_ac_e_org','EX_acald_e_org','EX_etoh_e_org','EX_lipa_cold_e_org',...
'EX_lipa_e_org','EX_colipa_e_org','EX_acolipa_e_org','EX_enlipa_e_org',...
'EX_kdo2lipid4_e_org', 'EX_hxa_e_org','EX_acser_e_org','EX_eca4colipa_e_org',...
'EX_his__L_e_org','EX_uri_e_org','EX_cytd_e_org','EX_xtsn_e_org','EX_ins_e_org1',...
'EX_adn_e_org','EX_thymd_e_org','EX_anhgm_e_org','EX_tyr__L_e_org',...
'EX_enter_e_org','EX_feenter_e_org','EX_phe__L_e_org','EX_pheme_e_org',...
'EX_ins_e_org'};
Com.ub(find(contains(Com.rxns,exchange_close))) = 0;
Com_save = Com; % this is the final matrix

% Calculate parsimonus maximal community growth rate
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
Com.lb(find(ismember(Com.rxns,'bio_sum'))) = sum(FBA(bio));
Com.c(:) = 0;
Com.c(find(contains(Com.rxns,'tr_s'))) = -1;
pFBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);

%% Add an equal biomass ratio if desired
Com = Com_save;
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(endsWith(Com.rxns,'_s'))) = 100;
Com.S(end+1:end+2,:) = 0;
Com.S(end-1,bio(1)) = -1;
Com.S(end-1,bio(2)) = 1;
Com.S(end,bio(2)) = -1;
Com.S(end,bio(3)) = 1;
Com.mets(end+1:end+2) = {'equal_growth'};
Com.b(end+1:end+2) = 0;
Com_equal = Com;

%% Phenotypic phase plane
Com = Com_save;
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
for i = 1:3
    Com.c(:) = 0;
    Com.c(bio(i)) = 1;
    FBA_self = cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    mu(:,i) = FBA_self(bio);
end
X_graph = meshgrid(0:mu(1,1)/19:mu(1,1),1:20);
Y_graph = [];
for i = 1:20
    Com = Com_save;
    Com.lb(find(ismember(Com.rxns,remove))) = 0;
    Com.ub(find(ismember(Com.rxns,remove))) = 0;
    Com.lb(bio(1)) = X_graph(1,i);
    Com.ub(bio(1)) = X_graph(1,i);
    Com.c(:) = 0;
    Com.c(bio(2)) = 1;
    FBAy_max=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    y_max(i) = FBAy_max(bio(2));
    Com.c(bio(2)) = -1;
    FBAy_min=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    y_min(i) = FBAy_min(bio(2));
    Y_graph(1:20,i) = linspace(y_min(i),y_max(i),20);
end

strain1 = [];
strain2 = [];
strain3 = [];

for i = 1:20
    for j = 1:20
        Com = Com_save;
        Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
        Com.lb(bio(1)) = X_graph(i,j);
        Com.lb(bio(2)) = Y_graph(i,j);
        Com.ub(bio(1)) = X_graph(i,j);
        Com.ub(bio(2)) = Y_graph(i,j);
        Com.c(:) = 0;
        Com.c(bio(3)) = 1;
        FBAppp=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if isempty(FBAppp) == 0
                strain1(i,j) = FBAppp(bio(1));
                strain2(i,j) = FBAppp(bio(2));
                strain3(i,j) = FBAppp(bio(3));
        else
            break
        end
    end
end
figure(1)


% colorbar
strain1 = strain1(:);
strain2 = strain2(:);
strain3 = strain3(:);
[k2,av1] = convhull(strain1,strain2,strain3,'Simplify',true);
trisurf(k2,strain1,strain2,strain3,'FaceColor',[0 0.4470 0.7410]);
xlabel('Strain 1 growth rate (h^{-1})')
ylabel('Strain 2 growth rate (h^{-1})')
zlabel('Strain 3 growth rate (h^{-1})')

%% Perform dFBA
met_list = upper(erase(erase(Com.rxns(find(endsWith(Com.rxns,'_s'))),'EX_'),'_s'));
% set-up initial parameters
    V=1;
    Yxs = 0.85;
    d = 0.1/3600;
    X0 = [0.1;0.1;0.1]; % initial biomass (g-DW/L)
    S0 = 10/180*1000; % inital glucose concentration (mmol/L)
    dt = 0.1; % space between iteration
    kd = 0.01 ;% decay
    km = 1; % m-m kinetics
    
% Scenario 1: Equal growth rate
    Com = Com_equal;
    i=1;
    S(i) = S0;
    X(:,i) = X0;
    while isempty(FBA)==0
        Glu_lb = round(-1*Com_save.lb(find(ismember(Com_save.rxns,'EX_glc__D_e_com')))*(S(i)/(S0+4/180)),3);
        Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))) = -1*Glu_lb;
        FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if sum(FBA(bio)) <= 0.0001
            break
        end
        X(:,i+1) = X(:,i).*exp((FBA(bio)-kd)*dt);
        S(i+1) = S(i)-Glu_lb*sum(X(:,i))*dt;
        i=i+1;
    end
    X_eq = X; 
    
% Scenario 2: Maximize community growth
    V=1;
    Yxs = 0.85;
    d = 0.1/3600;
    X0 = [0.1;0.1;0.1]; % initial biomass (g-DW/L)
    S0 = 10/180*1000; % inital glucose concentration (mmol/L)
    dt = 0.05; % space between iteration
    kd = 0.01 ;% decay
    km = 0.0001; % m-m kinetics
    Com = Com_save;
    Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
    Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
    % Reset variables and initialize
    X=0;
    S=0;
    i=1;
    S(i) = S0;
    X(1:3,i) = X0;
    met_conc(1:length(met_list),1) = 0;
    Glu_Vmax = Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))); % Depends on model
    % Begin dFBA
        while isempty(FBA)==0
            Glu_lb = round(-Glu_Vmax*(S(i)/(S(i)+4/180)),3); % 4 is from Bionumbers
            Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))) = -1*Glu_lb;
            FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
            if isempty(FBA)==1 | sum(FBA(bio)) <= 0.00001 | i > 480
                break
            end
            X(1:3,i+1) = X(:,i).*exp((FBA(bio)-kd)*dt);
            S(i+1) = max(0,S(i)-Glu_lb*sum(X(:,i))*dt);
            % Set-up organism ratio constraint, except for shared part,
            % growth and ATPM
            Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.lb(find(contains(Com.rxns,'_org3'))) = Com_save.lb(find(contains(Com.rxns,'_org3')))*X(3,i+1)/sum(X(:,i+1));            
            Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org3'))) = Com_save.ub(find(contains(Com.rxns,'_org3')))*X(3,i+1)/sum(X(:,i+1));            
            Com.ub(bio) = 100;
            Com.lb(bio) = 0;
            Com.ub(find(contains(Com.rxns,'ATPM'))) = 100;
            % shared is the minimal scenario captured by pFBA
            %Com.ub(find(contains(Com.rxns,'tr_s'))) = abs(pFBA(find(contains(Com.rxns,'tr_s'))));

            % Neglect the first hour, as we assume media is enough
            if i <= 10
                importation(1:length(met_list),1) = 0;
                importation(1:length(met_list),2) = 0; 
                importation(1:length(met_list),3) = 0;
            else
                importation(:,1) = abs(pFBA(find(contains(Com.rxns,'tr_s_org1_import'))));
                importation(:,2) = abs(pFBA(find(contains(Com.rxns,'tr_s_org2_import'))));
                importation(:,3) = abs(pFBA(find(contains(Com.rxns,'tr_s_org3_import'))));
            end   
            
            production(:,1) = abs(pFBA(find(contains(Com.rxns,'tr_s_org1_export'))));
            production(:,2) = abs(pFBA(find(contains(Com.rxns,'tr_s_org2_export'))));
            production(:,3) = abs(pFBA(find(contains(Com.rxns,'tr_s_org3_export'))));
            % Perform one iteration 
                for j = 1:length(met_list)
                    dx(j,i) = production(j,1)*X(1,i) + production(j,2)*X(2,i)+ production(j,3)*X(3,i)-...
                        importation(j,1)*X(1,i) - importation(j,2)*X(2,i)- importation(j,3)*X(3,i);
                    conc(j) = met_conc(j,i) + dx(j,i)*dt;
                    if conc(j) < 0
                       conc(j) = 0;
                    end
                    if i > 10
                        loc = find(contains(Com.rxns,strcat(met_list(j),'tr_s'))); % location of met
%                         Com.ub(loc(4)) = abs(pFBA(loc(4))*(conc(j)/(conc(j)+km))*X(1,i+1)/sum(X(:,i+1)));
%                         Com.ub(loc(5)) = abs(pFBA(loc(5))*(conc(j)/(conc(j)+km))*X(2,i+1)/sum(X(:,i+1)));
%                         Com.ub(loc(6)) = abs(pFBA(loc(6))*(conc(j)/(conc(j)+km))*X(3,i+1)/sum(X(:,i+1)));
                    end
                end
%             Com.ub(find(endsWith(Com.rxns,'_s'))) = 100;
            % re-introduce the kos
            Com.lb(find(ismember(mcs(:,MCS_chosen),-1)))=0;
            Com.ub(find(ismember(mcs(:,MCS_chosen),-1)))=0;
            met_conc(:,i+1) = conc;
            i=i+1;
        end    
    
