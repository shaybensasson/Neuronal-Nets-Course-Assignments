function nets = getNaturalImages(resizeN, numberOfNets)

% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
imagefiles = dir('C:\לימודים\תואר שני\ex2\images\cd01A\*.jpg');
nfiles = length(imagefiles);    % Number of files found

nets = cell(numberOfNets,1);
for i=1:numberOfNets
    currentfilename = imagefiles(i).name;
    currentimage = imread(currentfilename);
    
    
    currentimage = imresize(currentimage,[resizeN resizeN])
    currentimage = im2bw(currentimage,0.3);
    nets{i,1} = currentimage;
    %imagesc(BW);
    %colormap(gray)
end
end