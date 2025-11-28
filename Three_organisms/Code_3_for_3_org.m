clear;
cd 'C:\Users\alexa\OneDrive - University of Toronto\Documents\PhD_UofT\Paper\MCS_Com\Code\models'
addpath('C:\Users\alexa\OneDrive - University of Toronto\Documents\PhD_UofT\Paper\MCS_Com\Code\MCS');
addpath('C:\Users\alexa\OneDrive - University of Toronto\Documents\PhD_UofT\Paper\MCS_Com\Code\');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
% CNA path
addpath('C:\Users\alexa\CellNetAnalyzer')
load 'MCS_Com_iMLcore.mat'


%% Make the model
% To reduce the size of the problem, we will close as many EX rxns as
% possible. We will also close the bounds the tr_s to really ensure the
% metabolites wont be produced by the cell at all. in the desired strain

    EX_close = {'EX_ac_e_org','EX_acald_e_org','EX_etoh_e_org','EX_lipa_cold_e_org',...
        'EX_lipa_e_org','EX_colipa_e_org','EX_acolipa_e_org','EX_enlipa_e_org',...
        'EX_kdo2lipid4_e_org', 'EX_hxa_e_org','EX_acser_e_org','EX_eca4colipa_e_org',...
        'EX_his__L_e_org','EX_uri_e_org','EX_cytd_e_org','EX_xtsn_e_org','EX_ins_e_org1',...
        'EX_adn_e_org','EX_thymd_e_org','EX_anhgm_e_org','EX_tyr__L_e_org',...
        'EX_enter_e_org','EX_feenter_e_org','EX_phe__L_e_org','EX_pheme_e_org',...
        'EX_ins_e_org'};
    EX_shared = find(contains(Com.rxns,'EX_shared_'));
    Tr_s_org1 = find(contains(Com.rxns,strcat(upper(org1_dep),'tr_s_org1')));
    Tr_s_org2 = find(contains(Com.rxns,strcat(upper(org2_dep),'tr_s_org2')));
    Tr_s_org3 = find(contains(Com.rxns,strcat(upper(org3_dep),'tr_s_org3')));
    if ~exist('cnan','var')
        startcna(1);
        wait(parfevalOnAll(@startcna,0,1));
    end
    Com.ub(find(contains(Com.rxns,EX_close))) = 0;
    Com.lb(find(contains(Com.rxns,EX_close))) = 0;
    Com.lb(EX_shared) = 0;
    Com.ub(EX_shared) = 0;
    Com.ub(Tr_s_org1) = 0;
    Com.ub(Tr_s_org2) = 0;
    Com.ub(Tr_s_org3) = 0;
    
% Transform into a cnapy object
    cnap = CNAcobra2cna(Com);

% Add gprRules from Cobra to CNA
    cnap = CNAsetGenericReactionData_with_array(cnap,'geneProductAssociation',Com.grRules);

% Set each bounds to infinity, will not work otherwise
    max = max(Com.ub);
    cnap.reacMin(cnap.reacMin==-max) = -inf;
    cnap.reacMax(cnap.reacMax== max) =  inf;

% Find glucose and biomass
    idx.glc  = find(ismember(cnap.reacID,{'EX_glc__D_e_com'}));
    idx.bm  = find(ismember(Com.c,1)); 
% Find atpm to close them and prevent them from being considered
    idx.atpm_org1 = find(ismember(cnap.reacID,{'ATPM_org1'}));
    idx.atpm_org2 = find(ismember(cnap.reacID,{'ATPM_org2'}));
    idx.atpm_org3 = find(ismember(cnap.reacID,{'ATPM_org3'}));
    cnap.reacMin(idx.atpm_org1) = 0;
    cnap.reacMin(idx.atpm_org2) = 0;
    cnap.reacMin(idx.atpm_org3) = 0;
    
%% Uncomment to test feasibility 
    FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
    FBA(idx.bm)
    CNAplotPhasePlane(cnap,[idx.bm(1),find(ismember(Com.rxns,strcat('EX_',strcat(org1_dep,'_e_org1'))))])

%% Set up cMCS input
    modules{1}.sense = 'target';
    modules{1}.type = 'lin_constraints';
    modules{1}.V = zeros(1,cnap.numr);
    modules{1}.V(1,idx.Target(1)) = -1;
    modules{1}.v = -0.0001;

    modules{2}.sense = 'target';
    modules{2}.type = 'lin_constraints';
    modules{2}.V = zeros(1,cnap.numr);
    modules{2}.V(1,idx.Target(2)) = -1;
    modules{2}.v = -0.0001;

    modules{3}.sense = 'target';
    modules{3}.type = 'lin_constraints';
    modules{3}.V = zeros(1,cnap.numr);
    modules{3}.V(1,idx.Target(3)) = -1;
    modules{3}.v = -0.0001;  
    
    modules{4}.sense = 'desired';
    modules{4}.type = 'lin_constraints';
    modules{4}.V = zeros(3,cnap.numr);
    modules{4}.V(1,idx.bm(1)) = -1;
    modules{4}.V(2,idx.bm(2)) = -1;
    modules{4}.V(3,idx.bm(3)) = -1;
    modules{4}.v = zeros(3,1);
    modules{4}.v(:,1) = -0.01;
%% Specify is you want knock-ins (ki) and if some reactions should not be ko'ed
    rkoCost = ones(cnap.numr,1);
    Genes = ~cellfun(@isempty,Com.grRules);
    rkoCost(find(ismember(Genes,0))) = nan; % nan means its excluded from the possible solution
    rkoCost(find(contains(Com.rxns,'H2Ot_'))) = nan;
    rkoCost(find(contains(Com.rxns,'tr_s_'))) = nan;
    rkiCost = nan(cnap.numr,1); % Here no ki

maxSolutions = 1;
maxCost = 20; % Max number of kos
options.milp_solver = 'cplex';
options.milp_time_limit  = 120;
%options.milp_bigM = 10000; %Try true if you cannot find solutions
options.mcs_search_mode = 1;
verbose = 1;

%% Run cMCS and save outcome 

[mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
    CNAMCSEnumerator3(cnap,modules,rkoCost,rkiCost,maxSolutions,maxCost,options,verbose);

% [mcs, status, cmp_mcs, cmp_cnap, mcs_idx_cmp_full, obj] = ...
%     CNAgeneMCSEnumerator3(cnap,modules,{},{},maxSolutions,maxCost,gkoCost,gkiCost,{},options,verbose);


%% Saving the solutions and the model used
    Com.keep = Keep;
    save ('ComStrategy_3_Core.mat','mcs');
    save ('MCS_3Core_test.mat','Com')  
    mcs_test = mcs(:,1);
    cnap.reacMin(mcs_test == -1) = 0;
    cnap.reacMax(mcs_test == -1) = 0;
    %CNAplotPhasePlane(cnap,idx.bm)
    cnap.objFunc(idx.bm) = 1;

%% Compute the flux space to confirm MCSs
    FBA = CNAoptimizeFlux(cnap,[],[],2,0,0,[],[]); 
    FBA(idx.bm)