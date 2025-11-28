clear;
cd '/nfs/homes/tremb133/CommunityStuff/Ratio/model';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_bacc.mat
%remove = {'PPNDH_org1';'DRPA_org2';'PGK_org2';'TALA_org2';'TPI_org2'};
%remove = {'PPNDH_org1';'DRPA_org2';'FGLU_1_org2';'G6PDH2r_org2';'TALA_org2'};
%remove = {'ASPTA_org1';'TYRTA_org1';'ADSS_org2';'ALLTN_org2';'DRPA_org2';'TALA_org2'};
remove = {'ASPTA_org1';'ILETA_org1';'LEUTAi_org1';'PHETA1_org1';'TYRTA_org1';'VALTA_org1';'DRPA_org2';'PGK_org2';'TALA_org2';'TPI_org2'};
bio = find(Com.c==1); % locate objective function

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
Com.ub(find(contains(Com.rxns,'tr_s'))) = 100;
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
'EX_ins_e_org','EX_dha_e_org','EX_chor_e_org','EX_fol_e_org','EX_rib__D_e_org',...
'EX_6pgc_e_org','EX_ribflv_e_org','EX_3ump_e_org','EX_ump_e_org','EX_3cmp_e_org','EX_cmp_e_org',...
'EX_gsn_e_org','EX_3gmp_e_org','EX_gmp_e_org','EX_3amp_e_org','EX_amp_e_org',...
'EX_dtmp_e_org','EX_spmd_e_org'};
Com.ub(find(contains(Com.rxns,exchange_close))) = 0;
%Com.ub(find(endsWith(Com.rxns,'_s'))) = 100;
idx.starch = find(contains(Com.rxns,'EX_starch_e'));
idx.glc  = find(ismember(Com.rxns,{'EX_glc__D_e_com'}));
Com.lb(idx.starch) = -1;
Com.lb(idx.glc) = 0;

%% Adding dextrin uptake into E. coli
Com.S(end+1,:) = zeros(1,size(Com.S,2)); % add one mets
Com.S(:,end+1:end+2) = zeros(size(Com.S,1),2); % add two reactionsl
% locate metabolites
id.h2o_c = find(ismember(Com.mets,'h2o_c_org1'));
id.h_c = find(ismember(Com.mets,'h_c_org1'));
id.pi_c = find(ismember(Com.mets,'pi_c_org1'));
id.atp_c = find(ismember(Com.mets,'atp_c_org1'));
id.adp_c = find(ismember(Com.mets,'adp_c_org1'));
id.dextr_com = find(ismember(Com.mets,'dextrin_e_org2'));
Com.mets(id.dextr_com) = {'dextrin_e_com'};
id.dextr_c = size(Com.S,1);
id.glc_c = find(ismember(Com.mets,'glc__D_c_org1'));
% Maltodextrin glucosidase dextrin (MLTG6_org1)
Com.S(id.dextr_c,end) = -1;
Com.S(id.h2o_c,end) = -2;
Com.S(id.glc_c,end) = 2;
% Dextrin transport in via ABC transport (DEXTRINt2_org1)
Com.S(id.dextr_com,end-1) = -1;
Com.S(id.h2o_c,end-1) = -1;
Com.S(id.atp_c,end-1) = -1;
Com.S(id.adp_c,end-1) = 1;
Com.S(id.dextr_c,end-1) = 1;
Com.S(id.h_c,end-1) = 1;
Com.S(id.pi_c,end-1) = 1;

% Add reaction and metabolite names
Com.rxns(end+1) = {'DEXTRINt2_org1'};
Com.rxns(end+1) = {'MLTG6_org1'};
Com.rxnNames(end+1) = {'Dextrin transport in via ABC transporters'};
Com.rxnNames(end+1) = {'Maltodextrin glucosidase dextrin'};
Com.mets(end+1) = {'dextrin_c_org1'};
Com.metNames(end+1) = {'Dextrine'};
% Add lb, ub, c and b
Com.lb(end+1:end+2) = -100;
Com.ub(end+1:end+2) = 100;
Com.b(end+1) = 0;
Com.c(end+1:end+2) = 0;
Com_save = Com; % this is the final matrix


%% Code to look at all possible combinaison
% load Gene_mcs_iAF_iYO.mat
% for i=1:size(real_removal,2)
%     remove = real_removal(~cellfun(@isempty,real_removal(:,i)),i);
%     Com = Com_save;
%     Com.lb(find(ismember(Com.rxns,remove))) = 0;
%     Com.ub(find(ismember(Com.rxns,remove))) = 0;
%     FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
%     phen(:,i) = FBA(bio);
% end
%% Calculate parsimonus maximal community growth rate
Com.lb(find(ismember(Com.rxns,remove))) = 0;
Com.ub(find(ismember(Com.rxns,remove))) = 0;
FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
FBA(bio)
Com.lb(find(ismember(Com.rxns,'bio_sum'))) = sum(FBA(bio));
Com.c(:) = 0;
Com.c(find(contains(Com.rxns,'tr_s'))) = -1;
pFBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);

%% Phenotypic phase plane
Com = Com_save;
X_ppp(1) = 0.5;
    X_ppp(2) = 0.5;
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

%% dFBA for selfish and optimal
met_list = upper(erase(erase(Com.rxns(find(endsWith(Com.rxns,'_s'))),'EX_'),'_s'));
% optimal scenario
    % Here initially, for the first hour we will assume they have enough
    % things to import from the environment, then we calculate mu, allowing
    % us to calculate the X, which will be used to calculate the ub on imporation.
    % set-up initial parameters
    V=1;
    Yxs = 0.85;
    X0 = [0.1;0.1]; % initial biomass (g-DW/L)
    S0 =10/0.324; % inital starch concentration assuming 2 unit of glucose (mmol/L) 
    Dex0 = 0; % concentration of initial dextrose
    km_am = 1.2; %mmol/L
    Vmax_am = 0.8; % Yoneda et al,1975
    km_dex = 0.0014;
    Vmax_dex = 0.65;
    Vmax = 1000;
    dt = 0.5; % space between iteration
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
    Dex(i) = Dex0;
    X(1:2,i) = X0;
    amylase = 0;
    met_conc(1:length(met_list),1) = 0;
    alpha_am = find(ismember(Com.rxns,'AAMYL_1_org2'));
    dex_imp = find(contains(Com.rxns,'DEXTRINt2_org'));
    Com.lb(idx.starch) = -100; % keep it unlimiting until these is none
    % Begin dFBA
        while isempty(FBA)==0
            Com.c(:) = 0;
            Com.c(bio) = 1;
            Com.lb(find(ismember(Com.rxns,'bio_sum'))) = 0;
            aamyl_lb = amylase(i) * 100/1000/1000*60/0.324*S(i)/(S(i)+km_am); % rate of alpha amylase in mmol/L*h
            Com.ub(dex_imp) = Vmax_dex*Dex(i)/(Dex(i)+km_dex);
            FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
            if isempty(FBA)==1 | sum(FBA(bio)) <= 0.001 | Dex(i) == 0 | i > 1000 
                if i > 5 
                    break
                end
            end
            amylase(i+1) = amylase(i) + 6500/24*X(2,i)*dt;
            X(1:2,i+1) = X(:,i).*exp((FBA(bio)-kd)*dt);
            S(i+1) = S(i)-aamyl_lb*dt;
            if S(i+1) <= 0.0001
                S(i+1) = 0;
            end
            Dex(i+1) = Dex(i) + abs(S(i+1)-S(i))-sum(FBA(dex_imp))*sum(X(i))*dt;
            if Dex(i+1) <= 0.0001
                Dex(i+1) = 0;
            end
            % Set-up organism ratio constraint, except for shared part,
            % growth and ATPM
            Com.lb(find(contains(Com.rxns,'_org1'))) = Com_save.lb(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.lb(find(contains(Com.rxns,'_org2'))) = Com_save.lb(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org1'))) = Com_save.ub(find(contains(Com.rxns,'_org1')))*X(1,i+1)/sum(X(:,i+1));
            Com.ub(find(contains(Com.rxns,'_org2'))) = Com_save.ub(find(contains(Com.rxns,'_org2')))*X(2,i+1)/sum(X(:,i+1));
            Com.ub(bio) = 100;
            Com.lb(bio) = 0;
            Com.ub(find(contains(Com.rxns,'ATPM'))) = 400;
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
                    if i>=100
%                         loc_org1 = find(ismember(Com.rxns,strcat(met_list(j),'tr_s_org1_import')));
%                         loc_org2 = find(ismember(Com.rxns,strcat(met_list(j),'tr_s_org2_import')));
%                         Com.ub(loc_org1) = Vmax*(conc(j)/(conc(j)+km))*X(1,i+1)/sum(X(:,i+1));
%                         Com.ub(loc_org2) = Vmax*(conc(j)/(conc(j)+km))*X(2,i+1)/sum(X(:,i+1));
                        %Com.ub(find(endsWith(Com.rxns,'_s'))) = 100;
                        %Com.lb(find(endsWith(Com.rxns,'_s'))) = -100;
                    end
                end
            Com.lb(find(ismember(Com.rxns,remove))) = 0;
            Com.ub(find(ismember(Com.rxns,remove))) = 0;
            met_conc(:,i+1) = conc;
            %FBA(bio)
            i=i+1;
        end
    X_pres = X';
    S_pres = S'*0.324;
    Dex_pres = Dex'*0.324;