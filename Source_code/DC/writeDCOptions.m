function status = writeDCOptions(opt,inifile)
% write ini config file

ini = IniConfig();
ini.AddSections('Parameters');
ini.AddSections('General');
ini.AddSections('Hypothesis Space');
ini.AddSections('Initialization');
ini.AddSections('Detections');
ini.AddSections('Functions');
ini.AddSections('Discrete Optimization');
ini.AddSections('Continuous Optimization');
ini.AddSections('Miscellaneous');

opttmp=opt;

while ~isempty(fieldnames(opttmp))
    fnames=fieldnames(opttmp);
    fieldname=fnames{1};
    ini=parseField(ini,opt,fieldname);
    opttmp=rmfield(opttmp,fieldname);
end

status = ini.WriteFile(inifile);

end

function ini=parseField(ini,opt,fieldname)
    fvalue=getfield(opt,fieldname);
    if isstruct(fvalue)
        opttmp=fvalue;
        switch (fieldname)
            % Detections
            case 'detScale'
                while ~isempty(fieldnames(opttmp))
                    fnames=fieldnames(opttmp);
                    fieldname=fnames{1};
                    fvalue=getfield(opttmp,fieldname);
                    ini.AddKeys('Detections',fieldname,fvalue);
                    opttmp=rmfield(opttmp,fieldname);
                end                
            % Discrete Params
            case 'disOpt'
                while ~isempty(fieldnames(opttmp))
                    fnames=fieldnames(opttmp);
                    fieldname=fnames{1};
                    fvalue=getfield(opttmp,fieldname);
                    ini.AddKeys('Discrete Optimization',fieldname,fvalue);
                    opttmp=rmfield(opttmp,fieldname);
                end                
            % Continuous Params
            case 'conOpt'
                while ~isempty(fieldnames(opttmp))
                    fnames=fieldnames(opttmp);
                    fieldname=fnames{1};
                    if strncmp(fieldname,'enPar',5)
                        opttmp=rmfield(opttmp,fieldname);
                        continue;
                    end
                    fvalue=getfield(opttmp,fieldname);
                    ini.AddKeys('Continuous Optimization',fieldname,fvalue);
                    opttmp=rmfield(opttmp,fieldname);
                end                
        end
        % frames, special case
    elseif strcmpi(fieldname,'frames')
        ini.AddKeys('General','ff',fvalue(1));
        ini.AddKeys('General','lf',fvalue(end));
    elseif strcmpi(fieldname,'seqLength')
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
    case {'unaryFactor', ...
            'labelCost', ...
            'outlierCost', ...
            'persistenceFactor', ...
            'curvatureFactor', ...
            'slopeFactor', ...
            'proxcostFactor', ...
            'exclusionFactor', ...
            'pairwiseFactor' ...
          }
      %             'fidelityFactor', ...
      %             'segFactor', ...
      %             'goodnessFactor' ...
        sec='Parameters';
        
        % General
    case {'track3d',...
            'verbosity', ...	
            'mex', ...		
            'visOptim', ...	
            'met2d', ...		
            'keepHistory', ...	
            'cutToTA', ...	
            'randrun', ...	
            'remOcc', ...		
            'maxItCount', ...	
            'occ', ...		
            'minCPs', ...	
            'ncpsPerFrame', ...	
            'tau', ...		
            'borderMargin' ...            
         }
        sec='General';
    case {
            'nInitModels', ...
            'maxModels', ...       
            'nMaxAddExt', ...
            'nMaxAddMerged', ...
            'nAddRandomModels', ...
            'nAddModelsFromOutliers' ...           
            }
        sec='Hypothesis Space';
    case {
            'startFromEKF', ...
            'startFromPir', ...
            'startFromGT', ...
            'EKFDir', ...
            'DPDir' ...            
            }
        sec='Initialization';
    case {
            'detThreshold', ...
            'sigA', ...
            'sigB' ...
            }
        sec='Detections';
    case {
            'dataFunction', ... 
            'speedFunction', ...
            'fidFunction' ...
            }
        sec='Functions';
   
end
        

end