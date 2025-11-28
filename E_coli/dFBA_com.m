clear;
cd '/nfs/homes/tremb133/CommunityStuff/Ratio/model';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load MCS_Com_iAF
remove = {'CS_org1';'GND_org1';'HCYSMT_org1';'HCYSMT2_org1';'IPPS_org1';'MALS_org1';'TALA_org1';'TRDR_org1';...
    'CHORM_org2';'PPND_org2';'PPNDH_org2'};
%   remove = {'DRPA_org1','HCYSMT_org1','HCYSMT2_org1','PPNDH_org1','TALA_org1',...
%      'HCYSMT_org2', 'HCYSMT2_org2','PPM2_org2','PPND_org2','TALA_org2'};
bio = find(Com.c==1); % locate objective function
ratio = 0;
%% Calculate minimal parsimonious exchanges
% Split shared exchange and import to find lowest flux strategy
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
% Remove possibility to waste carbon just to get flux out
exchange_close = {'EX_ac_e_org','EX_acald_e_org','EX_etoh_e_org','EX_lipa_cold_e_org',...
'EX_lipa_e_org','EX_colipa_e_org','EX_acolipa_e_org','EX_enlipa_e_org',...
'EX_kdo2lipid4_e_org', 'EX_hxa_e_org','EX_acser_e_org','EX_eca4colipa_e_org',...
'EX_his__L_e_org','EX_uri_e_org','EX_cytd_e_org','EX_xtsn_e_org','EX_ins_e_org1',...
'EX_adn_e_org','EX_thymd_e_org','EX_anhgm_e_org','EX_tyr__L_e_org',...
'EX_enter_e_org','EX_feenter_e_org','EX_phe__L_e_org','EX_pheme_e_org',...
'EX_ins_e_org'};
Com.ub(find(contains(Com.rxns,exchange_close))) = 0;
%Com.ub(find(endsWith(Com.rxns,'_s'))) = 100;
Com_save = Com; % this is the final matrix

% Calculate parsimonus maximal community growth rate
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
Com.lb(find(ismember(Com.rxns,'bio_sum'))) = sum(FBA(bio));
Com.c(:) = 0;
Com.c(find(contains(Com.rxns,'tr_s'))) = -1;
pFBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);

%% dFBA for selfish and optimal
met_list = upper(erase(erase(Com.rxns(find(endsWith(Com.rxns,'_s'))),'EX_'),'_s'));
% optimal scenario
    % Here initially, for the first hour we will assume they have enough
    % things to import from the environment, then we calculate mu, allowing
    % us to calculate the X, which will be used to calculate the ub on imporation.
    % set-up initial parameters
    V=1;
    Yxs = 0.85;
    X0 = [0.1;0.4]; % initial biomass (g-DW/L)
    S0 = 10/180*1000; % inital glucose concentration (mmol/L)
    dt = 0.1; % space between iteration
    t0= 0;
    kd = 0.01 ;% decay
    km = 0.0001; % m-m kinetics
    Com = Com_save;
    Com.lb(find(ismember(Com.rxns,remove))) = 0;
    Com.ub(find(ismember(Com.rxns,remove))) = 0;
    % Reset variables and initialize
    X=0;
    S=0;
    i=1;
    S(i) = S0;
    X(1:2,i) = X0;
    met_conc(1:length(met_list),1) = 0;
    Glu_Vmax = Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))); % Depends on model
    % Begin dFBA
        while isempty(FBA)==0
            Com.c(:) = 0;
            Com.c(bio) = 1;
            Com.lb(find(ismember(Com.rxns,'bio_sum'))) = 0;
            Glu_lb = round(-Glu_Vmax*(S(i)/(S(i)+4/180/1000)),3); % 4 mg/L from monod
            Com.lb(find(ismember(Com.rxns,'EX_glc__D_e_com'))) = -1*Glu_lb;
            FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
            if isempty(FBA)==1 | sum(FBA(bio)) <= 0.0001 | i > 480 
                break
            end
            X(1:2,i+1) = X(:,i).*exp((FBA(bio)-kd)*dt);
            S(i+1) = S(i)-Glu_lb*sum(X(:,i))*dt;
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
            %Com.ub(find(contains(Com.rxns,'tr_s'))) = pFBA(find(contains(Com.rxns,'tr_s')));

            % Neglect the first 1 hour, as we assume media is enough
            if i <= 0
                importation(1:length(met_list),1) = 0;
                importation(1:length(met_list),2) = 0;            
            else
                importation(:,1) = pFBA(find(contains(Com.rxns,'tr_s_org1_import')));
                importation(:,2) = pFBA(find(contains(Com.rxns,'tr_s_org2_import')));
            end   
            
            production(:,1) = pFBA(find(contains(Com.rxns,'tr_s_org1_export')));
            production(:,2) = pFBA(find(contains(Com.rxns,'tr_s_org2_export')));
            % Perform one iteration 
                for j = 1:length(met_list)
                    dx(j,i) = production(j,1)*X(1,i) + production(j,2)*X(2,i)-...
                        importation(j,1)*X(1,i) - importation(j,2)*X(2,i);
                    conc(j) = met_conc(j,i) + dx(j,i)*dt;
                    if conc(j) < 0
                       conc(j) = 0;
                    end
                end
            Com.lb(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(ismember(Com.rxns,remove))) = 0;
            met_conc(:,i+1) = conc;
            i=i+1;
        end
%% Phenotypic phase plane
Com = Com_save;
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

%% Pheotypic phase plane for different ratios
% case 1 - 1:2 and 2:1
if ratio == 1
    Com = Com_save;
    X_ppp(1) = 10;
    X_ppp(2) = 0.1;
    Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
    Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
    Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
    Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
    Com.ub(bio) = 100;
    Com.lb(bio) = 0;
    Com.ub(find(contains(Com.rxns,'ATPM'))) = 100;
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
            Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
            Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
            Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
            Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
            Com.ub(bio) = 100;
            Com.lb(bio) = 0;
            Com.ub(find(contains(Com.rxns,'ATPM'))) = 100;
            Com.lb(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
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
            Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
            Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
            Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X_ppp(1)/sum(X_ppp);
            Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X_ppp(2)/sum(X_ppp);
            Com.ub(bio) = 100;
            Com.lb(bio) = 0;
            Com.ub(find(contains(Com.rxns,'ATPM'))) = 100;
            Com.lb(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(endsWith(Com.rxns,'_s'))) = 0;
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
end