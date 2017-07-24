function convvid2img(filename)
    workingDir = extractBefore(filename, '.');
    if(~exist(workingDir))
        mkdir(workingDir)
        mkdir(workingDir,'images')
    end
    reader = VideoReader(filename);
    ii = 1;

    while hasFrame(reader)
        img = readFrame(reader);
        filename = [sprintf('%04d',ii) '.jpg'];
        fullname = fullfile(workingDir,'images',filename);
        imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
        ii = ii+1;
    end

end