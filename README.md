# FunGreen-papersHooftman-
Codes belonging to two Manuscripts (papers hopefully)
This folder contains the codes used for two MS’s  with working title 
1.	Arrested Dispersal, Eutrophication and External Seed Pressure Severely Degrade Grassland Green Infrastructure.
2.	Higher Ecosystem Services Credits derive from targeted Grassland Green Infrastructure management.
Codes and MS’s are copyrighted to DAP Hooftman, Lactuca: Environmental Data Analyses and Modelling, Diemen, The Netherlands.

Data will be published elsewhere after publication.

This folder contains:
-	GIStatsMain: contains the main steering code and the numbering of the different traits (LisFunc)
-	MakeMeanTraits: collation of all segment and core variables, normalisation of all mean trait values against respective cores, including:
-	CollateTraits: the calculation of all segment mean traits from species composition bespoke where needed or via:
-	StandNum: standardised function for mean trait calculation per segment , where possibly used.
-	CalcStatsRadi: the actual statistics comparing cores with segments (“differences”), GLM’s of structural connectivity against mean trait values per segment (“Anova_Table”), and among scenario comparisons (“Distances_overview”). All per mean trait and with trait categories if applicable.
To be added:
-	PredictFitMain: main and collating module of predictions and sensitivities (ES paper) steering for x runs of:
-	PredictLoop: actual prediction and sensitivity calculations for ES Paper, including:
-	RegressFunc: the two tiered  GLM calculations used for calculation of marginal means

Danny Hooftman, 30-03-2020
