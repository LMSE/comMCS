clear;
cd '/nfs/homes/tremb133/CommunityStuff/Ratio/model';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load MCS_Com_iMLcore_OLD.mat
load ComStrategy_4_Core


%% List all reactions in the MCS and perform analysis
bio = find(ismember(Com.c,1));

for i = 1:size(mcs,2)
    Solutions(1:size(find(ismember(mcs(:,i),-1)),1),i) = Com.rxns(find(ismember(mcs(:,i),-1)));
    mcs_size(i) = length(find(~ismember(mcs(:,i),0)));
end
exchange_close = {'EX_ac_e_org','EX_acald_e_org','EX_etoh_e_org','EX_lipa_cold_e_org',...
'EX_lipa_e_org','EX_colipa_e_org','EX_acolipa_e_org','EX_enlipa_e_org',...
'EX_kdo2lipid4_e_org', 'EX_hxa_e_org','EX_acser_e_org','EX_eca4colipa_e_org',...
'EX_his__L_e_org','EX_uri_e_org','EX_cytd_e_org','EX_xtsn_e_org','EX_ins_e_org1',...
'EX_adn_e_org','EX_thymd_e_org','EX_anhgm_e_org','EX_tyr__L_e_org',...
'EX_enter_e_org','EX_feenter_e_org','EX_phe__L_e_org','EX_pheme_e_org','EX_ins_e'};
Com.lb(find(contains(Com.rxns,exchange_close))) = 0;
Com.ub(find(contains(Com.rxns,exchange_close))) = 0;
Com_save = Com;

for i = 1:size(mcs,2)
    Com = Com_save;
    Com.lb(find(ismember(mcs(:,i),-1)))=0;
    Com.ub(find(ismember(mcs(:,i),-1)))=0;
    Com.lb(find(contains(Com.rxns,'EX_shared'))) = 0;
    Com.ub(find(contains(Com.rxns,'EX_shared'))) = 0;
    Com.c(:) = 0;
    Com.c(bio(1)) = 1;
    FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
    if sum(FBA(bio)~=0) == 3
        Com.c(:) = 0;
        Com.c(bio(2)) = 1;
        FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
        if sum(FBA(bio)~=0) == 3
            Com.c(:) = 0;
            Com.c(bio(3)) = 1;
            FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
            if sum(FBA(bio)~=0) == 3
                works(i) = 1;
                bio_val(1:3,i) = FBA(bio);
            else
                works(i) = 0;
                bio_val(1:3,i) = FBA(bio);
            end
        else
            works(i) = 0;
            bio_val(1:3,i) = FBA(bio);      
        end
    else 
        works(i) = 0;
        bio_val(1:3,i) = FBA(bio);
    end
end