#!/bin/bash

# compute substitution models of different partitions


# -AIC -BIC : find AIC, BIC
# -tr 8: use 8 threads
# -a : calculate a model averaged phylogeny
# -f : Include models with unequals base frecuencies.
# -g numberOfRateCategories Include models with rate variation among sites and sets the number of categories. Usually 4 categories are enough.
# -i : Include models with a proportion invariable sites.

# model for variant sites dataset
java -jar ~/jmodeltest-2.1.10/jModelTest.jar -AICc -tr 4 -g 4 -i -f  -a -d island_mainland_all_trimmed.fasta -o data/output/model_selection/all_SNP &

# model for full  dataset
java -jar ~/jmodeltest-2.1.10/jModelTest.jar -AICc -tr 4 -g 4 -i -f  -a -d island_mainland_all.fasta -o data/output/model_selection/full_alignment &

# model for variant sites with no missing data
java -jar ~/jmodeltest-2.1.10/jModelTest.jar -AICc -tr 4 -g 4 -i -f  -a -d island_mainland_all_trimmed_267.fasta -o data/output/model_selection/267 &

# model for variant and invariant sites with no missing data
java -jar ~/jmodeltest-2.1.10/jModelTest.jar -AICc -tr 4 -g 4 -i -f  -a -d island_mainland_all_alignment_complete_deletion.fasta -o data/output/model_selection/complete_deletion &



# compute ML phylogenetic tree

# save as a phylip file in aliview and then:

# do 1000 bootstraps of a kimura model 
#"New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0."
#Guindon S., Dufayard J.F., Lefort V., Anisimova M., Hordijk W., Gascuel O.
#Systematic Biology, 59(3):307-21, 2010.
phyml -i island_mainland_all_trimmed_267.phy -q -b 1000 -m K80 &




# -b: 1000 bootstraps
# -m: use GTR model estimated with jmodeltest
# -f: the character frequencies are determined with ML
# -v: proportion of invariable sites estimated with ML
# -c: number of substitution rate categories
# -a: estimate gamma parameter with ML
# -o tlr : tree topology (t), branch length (l) and rate parameters (r) are optimised.
# -s: best of NNI and SPR search. SPR = subtree pruning and regrafting (SPR) = generally better
phyml -i island_mainland_all_alignment_complete_deletion.phy -b 1000 --run_id GTR+I+G -m GTR -f m -v e -c 4 -a e --no_memory_check -o tlr -s BEST &
