function outpts=selectPointsSubset(alldpoints,selection)
% TODO



    Field = fieldnames(alldpoints);

        for iField = 1:length(Field)
            fcontent=alldpoints.(char(Field(iField)));
            outpts.(char(Field(iField)))=fcontent(selection);
        end

end