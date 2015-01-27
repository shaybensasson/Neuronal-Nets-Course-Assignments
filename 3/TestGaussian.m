figure;

hSize = [10 10];
SIGMAS = [1; 5];
subplot(length(SIGMAS)+1,1,1)
%%# Read an image
I = imread('peppers.png');
imshow(I)

for i=1:length(SIGMAS)
    %# Create the gaussian filter with hsize and sigma
    
    G = fspecial('gaussian',hSize,SIGMAS(i));
    %# Filter it
    Ig = imfilter(I,G,'same');
    %# Display
    subplot(length(SIGMAS)+1,1,i+1)
    imshow(Ig)
end

