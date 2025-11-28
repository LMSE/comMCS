clear;
cd '/nfs/homes/tremb133/CommunityStuff/ComMCS';
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio');
addpath('/nfs/homes/tremb133/CommunityStuff/Ratio/model');
% Cplex path
addpath('/nfs/homes/tremb133/CPLEX_Studio128/cplex/matlab/x86-64_linux') 
load /nfs/homes/tremb133/CommunityStuff/Ratio/model/MCS_Com_iAF.mat
load October_iAF_MCS  % your MCS file

bio = find(ismember(Com.c,1));
Com_save = Com;

%% Find what the rxns in the MCS are
for i = 1:size(mcs,2)
    Solutions(1:size(find(ismember(mcs(:,i),-1)),1),i) = Com.rxns(find(ismember(mcs(:,i),-1)));
end

%% Test if the genes work 
% Here we just need to make sure that only the genes of one strains are
% removed
  %for i = 3313
  real_removal = {};
  mcs_val = [];
  for i=1:size(mcs,2)
  %for i = 2157:2172
     %  Find the genes associated to each rxns
     rxnKOs = find(ismember(mcs(:,i),-1));
     if size(rxnKOs,1) < 10
        geneKOs = Com.grRules(find(ismember(mcs(:,i),-1)));
        for m = 1:size(rxnKOs,1)
            if contains(Com.rxns(rxnKOs(m)),"_org1")
                organism(m) = 1;
            elseif contains(Com.rxns(rxnKOs(m)),"_org2")
                organism(m) = 2;
            else
                organism(m) = 0;
            end
        end
        % From these get list of possibilities
        kill_genes = ["temp"];
        gene_org = 0;
        for j = 1:size(geneKOs,1)
            size_genes = size(kill_genes,2);
            if contains(geneKOs(j),'or')
                to_remove = split(geneKOs(j)," or ")';
                for k = 1:length(to_remove)
                    kill_genes(:,end+1) = to_remove(k);   
                end
            elseif contains(geneKOs(j),'and')
                possible = split(geneKOs(j)," and ");
                kill_genes(:,end+1) = possible(1);
                current_size = size(kill_genes,1);
                for l = 2:length(possible)
                    kill_genes(end+1:end+current_size,:) = kill_genes(1:current_size,:);
                    kill_genes(current_size+l-1:size(kill_genes,1),end) = possible(l);
                end
            else
                kill_genes(:,end+1) = geneKOs(j);
            end
            gene_org(1,size_genes+1:size(kill_genes,2)) = organism(j);
        end
        kill_genes(:,1) = [];
        gene_org(:,1) = [];
        
 % Find rxns that will be ko'ed with each possible solutions
   clear j k l m rmvRxns
   
   rmvRxns(1,1:size(kill_genes,1)+1) = 0;
        for j=1:length(Com.grRules)
            if contains(Com.rxns(j),"_org1")
                type = 1;
            elseif contains(Com.rxns(j),"_org2")
                type = 2;
            else
                type = 0;               
            end
            rules = Com.grRules(j);
            rmvRxns(end+1,:) = 0;
            rmvRxns(end,1)=j;
            if type == 0
                continue
            end
            if cellfun(@isempty,Com.grRules(j)) == 1
                continue
            elseif contains(rules,'or')
                all_genes = split(rules," or ")';
                for k = 1:size(kill_genes,1)
                    if length(find(ismember(gene_org(find(ismember(kill_genes(k,:), all_genes))),type))) == length(all_genes)
                       rmvRxns(end,k+1) = 1;
                    end
                end
            elseif contains(rules,'and')
                all_genes = split(rules," and ")';
                for k = 1:size(kill_genes,1)
                    if length(find(ismember(gene_org(find(ismember(kill_genes(k,:), all_genes))),type))) >= 1
                       rmvRxns(end,k+1) = 1;
                    end
                end
            else
                for k = 1:size(kill_genes,1)
                    unique_kill = unique(kill_genes(k,:));
                    if length(find(ismember(gene_org(find(ismember(kill_genes(k,:), rules))),type))) == 1 
                        rmvRxns(end,k+1) = 1;
                    end
                end
            end
        end
     rmvRxns(1,:) = [];
     % Finding if the ko's kill the cell
     clear j k works
     for j = 2:size(rmvRxns,2);
         Com = Com_save;
         Com.lb(find(ismember(rmvRxns(:,j),1))) = 0;
         Com.ub(find(ismember(rmvRxns(:,j),1))) = 0;
         FBA=cplexlp(-Com.c,[],[],Com.S,Com.b,Com.lb,Com.ub);
         works(j-1) = (sum(FBA(bio)) ~= 0);
     end
     MCS_valid(1,i) = sum(works);
     MCS_valid(2,i) = length(unique(kill_genes));
     if sum(works) > 0
        for j = 1:size(kill_genes,1)
            real_removal(1:size(find(ismember(rmvRxns(:,j+1),1)),1),end+1) = Com.rxns(find(ismember(rmvRxns(:,j+1),1)));
            mcs_val(end+1) = i;
        end
     end
     end
  end
  save ('Gene_mcs_iAF_keio_Oct.mat','real_removal');
  save ('MCS_valid_iAF_keio_Oct.mat','MCS_valid');
  save ('mcs_val_Oct.mat','mcs_val');