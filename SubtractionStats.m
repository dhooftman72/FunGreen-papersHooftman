
List = {'LogNectar','PollinatorVisits_cover','CarbonLeaf','Release_height_Cover','DryWeightLeaf',...
    'CarbonRoot','DryWeightRoot','RootingDepth','LogSpecies','SpeciesShannon','ColourShannon',...
    'Life_form_Diversity','Landscape_Diversity'}

for i = 1:13
Line1(i,1:3) = Predictions.(genvarname([char(List(i))])).Single_trait.Correlations(1,2:4);
Line1(i,4:8) = Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,1:5);
Line2(i,1:3) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,4:6);
Line2(i,4) = Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,9);
Line2(i,5:7) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,13:15);
Line2(i,8) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,13);
Line2(i,9:11) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,22:24);
Line2(i,12) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,17);
Line2(i,13:15) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,31:33);
Line2(i,16) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,21);
Line2(i,17:19) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,40:42);
Line2(i,20) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,25);

Normal(i,1) = Output_Anova.(genvarname([char(List(i))])).Single_trait.MeanCattle.Anova_Table(15,6);
end

List = {'LogSpecies','SpeciesShannon'};
for i = 1:2
Line1Single(i,1:3) = Predictions.(genvarname([char(List(i))])).Single_trait.Correlations(1,2:4);
Line1Single(i,4:8) = Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,1:5);
Line2Single(i,1:3) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,4:6);
Line2Single(i,4) = Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,9);
Line2Single(i,5:7) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,13:15);
Line2Single(i,8) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,13);
Line2Single(i,9:11) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,22:24);
Line2Single(i,12) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,17);
Line2Single(i,13:15) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,31:33);
Line2Single(i,16) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,21);
Line2Single(i,17:19) = Predictions.(genvarname([char(List(i))])).Sensitivities.Countries(1,40:42);
Line2Single(i,20) =Predictions.(genvarname([char(List(i))])).Sensitivities.All(1,25);

NormalSingle(i,1) = Output_Anova.(genvarname([char(List(i))])).Single_trait.MeanCattle.Anova_Table(16,6);
end

