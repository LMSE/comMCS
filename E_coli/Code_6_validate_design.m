clear;
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load To_test_iAF_keio_Oct.mat
load MCS_Com_iAF
load October_iAF_MCS.mat 
load Gene_mcs_iAF_keio_Oct.mat
load MCS_valid_iAF_keio_Oct.mat
load mcs_val_Oct.mat
bio = find(ismember(Com.c,1));
Com_save = Com;

%% Find what the rxns in the MCS are
for i = 1:length(to_test)
    Solutions(1:length(real_removal(find(~cellfun(@isempty,real_removal(:,to_test(i)))),to_test(i))),i) = real_removal(find(~cellfun(@isempty,real_removal(:,to_test(i)))),to_test(i));
end


%% Test if the genes work 
% Here we just need to make sure that only the genes of one strains are
% removed
clear i
  %for i = 1:30
  for i=1:size(to_test,2)
     %  Find the genes associated to each rxns
     rxnKOs = find(mcs(:,mcs_val(to_test(i)))==-1);
     geneKOs = Com.grRules(rxnKOs);
     kill_genes = ["temp"];
     gene_org = 0;
     flag = 0;
        for j = 1:size(geneKOs,1)
            size_genes = size(kill_genes,2);
            if contains(geneKOs(j),'))') == 1
               flag = 1;
            end
        end
     if flag == 1
        good(i) = 0;
     else
        good(i) = 1;
     end
  end
final_candidate = [];
for i=1:length(good)
    if good(i) == 1
        final_candidate(end+1) = to_test(i);
    end
end
cd '/nfs/homes/tremb133/CommunityStuff/ComMCS';
save('final_candidate_october_2024.mat','final_candidate')