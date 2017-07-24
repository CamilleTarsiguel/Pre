function homedir=getHomeFolder()
% home directory

    homedir='/home/amilan';
    if ispc
        homedir='D:';
    end
    if exist('/home/h3','dir')        
        homedir='/home/h3/amilan';
    end
end