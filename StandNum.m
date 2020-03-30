function [Selection,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi)
Selection.maxtraits(Trans.maxi) = Trans.Traitmax; 
name1 = ['Input.',Trans.naming,'(:,',num2str(Trans.tf(1)),')'];
RichnessTrait = zeros(Trans.nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(Trans.nrspeciesSegment,Trans.Traitmax);
StrengthTrait = ones(Trans.nrspeciesSegment,Trans.Traitmax);
if Trans.Traitmax == 1
    for x = 1:1:Trans.nrspeciesSegment
        PresList = find(eval(name1)== SpeciesSegment(x,1));
        name2 = ['Input.',Trans.naming,'(PresList,',num2str(Trans.tf(2)),')'];
        name3 = ['Input.',Trans.naming,'(PresList,',num2str(Trans.tf(3)),')'];
        if isempty(PresList) ~= 1
            RichnessTrait_T = (eval(name2));
            CoverTrait_T = SpeciesSegment(x,2);
            if isnan(Trans.tf(3))~= 1
                Strenght = (eval(name3));
            else
                Strenght =  1;
            end 
            RichnessTrait(x) = ((sum((bsxfun(@times, (RichnessTrait_T),Strenght))))./sum(Strenght));
            CoverTrait(x) = ((sum((bsxfun(@times, (CoverTrait_T),Strenght))))./sum(Strenght));
        else
            RichnessTrait(x) = 9999;
            CoverTrait(x) = 9999;
            StrengthTrait(x)= 9999;
        end
        clear RichnessTrait_T CoverTrait_T
    end
    clear PresList
    clear x
    [RichnessTrait,CoverTrait,~] = CleanTrait(RichnessTrait,CoverTrait,StrengthTrait);
    to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
    RichnessTrait =  RichnessTrait(to_take);
    CoverTrait =  CoverTrait(to_take);
    Selection.nrspecs = length(RichnessTrait);
    Selection.Total_cover = sum(CoverTrait); %
    Selection.CoverTrait = CoverTrait; %   
    RichnessTrait  = reshape(RichnessTrait ,length(RichnessTrait),1);
    CoverTrait = reshape(CoverTrait,length(CoverTrait),1);
    Nr_covers =  (bsxfun(@times, RichnessTrait, CoverTrait));
    
    if strcmp(Trans.which,'Both') || strcmp(Trans.which,'Prop')
        SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(Trans.segmentnr,1)  = mean(RichnessTrait);
        SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(Trans.segmentnr,1) =  length(RichnessTrait);
    end
    if strcmp(Trans.which,'Both') || strcmp(Trans.which,'Cover')
        SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(Trans.segmentnr,1)  =  sum(Nr_covers)./sum(CoverTrait);
        SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(Trans.segmentnr,1) = length(RichnessTrait);
    end
elseif Trans.Traitmax > 1 % categorical
    count = 0;
    for traittype =  Trans.TraitsToDo
        %Create Yes/no matrices per Life form
        count = count + 1;
        for x = 1:1:Trans.nrspeciesSegment
            PresList = find(eval(name1)== SpeciesSegment(x,1));
            name2 = ['Input.',Trans.naming,'(PresList,',num2str(Trans.tf(2)),')'];
            name3 = ['Input.',Trans.naming,'(PresList,',num2str(Trans.tf(3)),')'];
            if isempty(PresList) ~= 1
                Form = (eval(name2));
                 if isnan(Trans.tf(3))~= 1
                    Strenght = (eval(name3));
                else
                    Strenght = 1;
                end
                Pres = find(Form == traittype);
                if isempty(Pres) ~= 1
                    RichnessTrait(x,count) = 1;
                    CoverTrait(x,count) = SpeciesSegment(x,2);
                    StrengthTrait(x,count)= mean(Strenght);
                end
            else
                RichnessTrait(x,count) = 9999;
                CoverTrait(x,count)= 9999;
                StrengthTrait(x,count)= 9999;
            end
        end
        clear PresList
        clear Pres
        clear Forms
        clear x
    end
    [RichnessTrait,CoverTrait,StrengthTrait] = CleanTrait(RichnessTrait,CoverTrait,StrengthTrait);
    % Correct for strength of trait determination
    RichnessTrait = RichnessTrait.*StrengthTrait;
    CoverTrait = CoverTrait.*StrengthTrait;
    % calculate averages
    
    if strcmp(Trans.which,'Both') || strcmp(Trans.which,'Prop')
        SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(Trans.segmentnr,:) = sum(RichnessTrait,1);
        SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(Trans.segmentnr,:) =  SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(Trans.segmentnr,:);
    end
    if strcmp(Trans.which,'Both') || strcmp(Trans.which,'Cover')
        SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(Trans.segmentnr,:)  = sum(CoverTrait,1);
        SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(Trans.segmentnr,:) =  SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(Trans.segmentnr,:);
    end
end
end

function [RichnessTrait,CoverTrait,StrengthTrait] = CleanTrait(RichnessTrait,CoverTrait,StrengthTrait)
 % Remove species with no data
    RemoveDisp = find(RichnessTrait(:,1) == 9999);
    RichnessTrait(RemoveDisp,:) = [];
    CoverTrait(RemoveDisp,:) = [];
    StrengthTrait (RemoveDisp,:) = [];

    %check for zero entries
    listzeros = find(mean(RichnessTrait,2) == 0);
    RichnessTrait(listzeros,:) = [];
    CoverTrait(listzeros,:)= [];
    StrengthTrait(listzeros,:)= [];
end