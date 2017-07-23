function opt=parseOptimOpts(opt,optfile)


    xDoc=xmlread(optfile);
    
    % DISCRETE
    disOptTag=xDoc.getElementsByTagName('disOpt');
    
    
    da = char(disOptTag.item(0).getAttribute('alg'));
    opt.disOpt.alg = getDisAlg(da);
    infParamsTag = disOptTag.item(0).getElementsByTagName('infParams');
    opt.disOpt.initSimple = str2double(infParamsTag.item(0).getAttribute('initSimple'));
    
    theAttributes = getAttributes(infParamsTag.item(0));
    numPars = getLength(theAttributes);
    for n=1:numPars
        %Suggestion of Adrian Wanner
        str = toCharArray(toString(item(theAttributes,n-1)))';
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        p=sscanf(attr_name,'p%d');

        if ~isempty(p)
            opt.disOpt.infParam(p)=str2double(infParamsTag.item(0).getAttribute(attr_name));
        end
    end

    % CONTINUOUS
    conOptTag=xDoc.getElementsByTagName('conOpt');
    
    
    opt.conOpt.alg = char(conOptTag.item(0).getAttribute('alg'));
    infParamsTag = conOptTag.item(0).getElementsByTagName('infParams');
%     opt.disOpt.initSimple = str2double(infParamsTag.item(0).getAttribute('initSimple'));
%     
    theAttributes = getAttributes(infParamsTag.item(0));
    numPars = getLength(theAttributes);
    for n=1:numPars
        %Suggestion of Adrian Wanner
        str = toCharArray(toString(item(theAttributes,n-1)))';
        k = strfind(str,'=');
        attr_name = str(1:(k(1)-1));
        evalstr= ...
            sprintf('opt.conOpt.%s=str2double(infParamsTag.item(0).getAttribute(''%s''));',attr_name,attr_name);
        eval(evalstr);
    end
    
end

function dalg = getDisAlg(da)
switch (da)
    case 'TRWS'
        dalg=1;
    case 'QPBO'
        dalg=2;
    case 'MQPBO'
        dalg=3;
    case 'ICM'
        dalg=4;
    case 'TRWSi'
        dalg=5;
    case 'FastPD'
        dalg=6;
    case 'BP'
        dalg=7;
    case 'TRBP'
        dalg=8;
    case 'LF'
        dalg=9;


end
end
