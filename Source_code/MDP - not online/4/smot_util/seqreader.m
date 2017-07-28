classdef seqreader < handle
    properties(SetAccess = 'private', Hidden = true)
        seqPath;
        seqExt;
        seqName;
        seqNameRaw;
        imageStartId;
        imageCurrentId;
        idStart;
        idNumLen;
        
    end
    
    properties(SetAccess = 'private', GetAccess = 'public')
        Width;
        Height;
        NumOfFrames;
        CurrentPos;
        CurrentImageName;
    end
        
    methods
        % constructor
        function this = seqreader(inpath)
            
            % find if it is a folder or file
            ftype = exist(inpath,'file');
            
            if ftype==2         % file
                [this.seqPath,this.seqName,this.seqExt] = fileparts(inpath);
                
            elseif ftype==7     % folder
                
                % load file based on the operating system
                % thanks to Sebastiano Vascon
                if ispc==1
                    res2 = [ dir(fullfile(inpath, '*.jpg')) ; dir(fullfile(inpath, '*.png')) ;
                        dir(fullfile(inpath, '*.jpeg')) ; dir(fullfile(inpath, '*.ppm')) ];
                    res=[];
                    
                    for i=1:size(res2,1)
                        res=[res; res2(i).name];
                    end
                    
                    this.seqPath = inpath;
                    this.seqName = res(1,1:end-4);
                    this.seqExt = res(1,end-3:end);
                else
                    
                    [s,res] = system(['ls -1 ' inpath ' |egrep "\.png$|\.jpg$|\.jpeg$|\.ppm$"']);%
                    % get the first file
                    ind = strfind(res,char(10));
                    
                    this.seqPath = inpath;
                    this.seqName = res(1:ind(1)-5);
                    this.seqExt = res(ind(1)-4:ind(1)-1);
                    35;
                end
            else
                error('Path does not exist!');               
            end
            
            % find the sequence part of the file            
            this.idStart = regexp(this.seqName, '[0-9]+$');
            this.idNumLen = length(this.seqName)-this.idStart+1;
            this.seqNameRaw = this.seqName(1:this.idStart-1);
            this.imageStartId = str2double(this.seqName(this.idStart:end));
            this.imageCurrentId = this.imageStartId;
            
            
            % generate the file name
            imageName = [this.seqPath '/' this.seqNameRaw num2str(this.imageStartId,['%0' int2str(this.idNumLen) 'd']) this.seqExt];
                        
            % read the image
            frame = imread(imageName);
            
            [this.Height this.Width D] = size(frame);
            
            this.CurrentImageName = imageName;
            this.CurrentPos = 0;
            
            % find number of images
            imageList = dir([this.seqPath '/' this.seqNameRaw  '*' this.seqExt]);
            this.NumOfFrames = size(imageList,1);
            35;
        end
        
        % get next frame
        function frame = grabFrame(this)
            imageName = [this.seqPath '/' this.seqNameRaw num2str(this.imageCurrentId,['%0' int2str(this.idNumLen) 'd']) this.seqExt];
            
            
            % read the image
            frame = imread(imageName);
            
            this.CurrentImageName = imageName;
            this.CurrentPos = this.CurrentPos + 1;
            this.imageCurrentId = this.imageCurrentId + 1;            
            
      
        end        
%         % seeks to the given frame. not accurate
        function seek(this,pos)
            this.CurrentPos = pos;
            this.imageCurrentId = this.imageStartId + pos;            
        end
        
%         % needs destructor !!!
%         function delete(this)
%         end
    end
end