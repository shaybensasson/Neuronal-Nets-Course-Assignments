function patterns = getNaturalImages(NodesPerRowCol, totalPatterns, BW_CONVERTION_BINARY_TH)

% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
folder = '.\images\cd01A\';

imagefiles = dir([folder '*.jpg']);
nfiles = length(imagefiles);    % Number of files found
if totalPatterns<=nfiles
    patterns = cell(totalPatterns,1);
    for i=1:totalPatterns
        currentfilename = imagefiles(i).name;
        
        currentimage = imread([folder currentfilename]);
        
        
        currentimage = imresize(currentimage,[NodesPerRowCol NodesPerRowCol]);
        currentimage = im2bw(currentimage,BW_CONVERTION_BINARY_TH);
        
        %imagesc(currentimage);
        
        currentimage= double(currentimage);
        currentimage(currentimage==0)=-1;
      
        patterns{i} = currentimage;
                
        %colormap(gray)
    end
end
end