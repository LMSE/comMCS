# Engineering interdependence in synthetic microbial communities
In this repository, you will find every code and ressources used to reproduced the data presented in our paper.\
Although the codes presented are in Matlab, the result can be reproduced with CNApy.

## Folder description
**E_coli** Data generated in the wet lab\
**Three_organisms** Results of each of the analysis\
**Bacc_coli** These codes should allow you to recreate your own thermodynamic models which are going to be the input of all further analysis. We used iML1515 directly retrieved from the BiGG database\
**Narringenin** Contains the code used to directly take the data of what is feasible or not (Figure 3) and transforms it in way that makes it amenable for Escher maps.
**Models** Contains all the models used to generate the community models\


## Codes description
All foldiers contain the same series of code doing the same task\
**Code_1** Create a community model\
**Code_2** Adds the shared compartiment allowing for the exchange of metabolites\
**Code_3** Sets up and run the MCS algorithm\
**Code_4** Converts the MCS solution into genes and validate that the community can still grow\
**Code_5** Test solutions on individual knock out to ensure no growth is possible with individual models \
**Code_6** Allows for the analysis of all valid solution
**dFBA** Code used to perform phenotypic phase plane analysis and dFBA
Question? alexandre.tremblay@mail.utoronto.ca









