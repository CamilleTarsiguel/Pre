function status = writeSceneOptions(sceneInfo,inifile)
% write ini config file

ini = IniConfig();
ini.AddSections('Source');
ini.AddSections('Miscellaneous');

opttmp=sceneInfo;

while ~isempty(fieldnames(opttmp))
    fnames=fieldnames(opttmp);
    fieldname=fnames{1};
    ini=parseField(ini,sceneInfo,fieldname);
    opttmp=rmfield(opttmp,fieldname);
end

status = ini.WriteFile(inifile);

end

function ini=parseField(ini,sceneInfo,fieldname)
fvalue=getfield(sceneInfo,fieldname);
if strcmpi(fieldname,'frameNums')
    % ignore
else
    sec=getSectionName(fieldname);
    ini.AddKeys(sec,fieldname,fvalue);
end
end

function sec=getSectionName(fieldname)
sec='Miscellaneous';
switch(fieldname)
    % Parameters (optimized)
    case {'dataset', ...
            'sequence', ...
            'imgFolder', ...
            'frameRate', ...
            'detfile', ...
            }
        sec='Source';
end


end