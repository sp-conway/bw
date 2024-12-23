function morphedImage = RDCL_MorphImages(morphProp, image1, image2, varargin)

%
% Function: RDCL_MorphImages
%
% Morphs 2 or more grayscale images. All images must be same size.
%   Simple linear combination of pixels.
%
% Arguments:
%   morphProp = 1xn vector of proportion of each dot pattern present in
%     morph. If doesn't sum to 1, will be rescaled. n = number of images.
%   image1 = the first image
%   image2 = the second image.
%   varargin:
%     imageN = more images.
%
% Return:
%   morphedImage = the morphed image.
%
% Example:
%   newImage = RDCL_MorphImages([.5 .5], image1, image2);
%
% Notes:
%
% Andrew Cohen
% 6/26/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Make sure proportions sum to 1
morphProp = abs(morphProp);
morphProp = morphProp./sum(morphProp);

% Error checking
gridSize_x = size(image1, 2);
gridSize_y = size(image1, 1);

% image1
if size(image1, 3) ~= 1
    error('RDCL_MorphImages all images must be grayscale.');
end

% image2
if gridSize_x ~= size(image2, 2)
    error('RDCL_MorphImages all grid sizes must be equal.');
end

if gridSize_y ~= size(image2, 1)
    error('RDCL_MorphImages all grid sizes must be equal.');
end

if size(image2, 3) ~= 1
    error('RDCL_MorphImages all images must be grayscale.');
end

% Other images
for k = 1:1:length(varargin)
    if gridSize_x ~= size(varargin{k}, 2)
        error('RDCL_MorphImages all grid sizes must be equal.');
    end

    if gridSize_y ~= size(varargin{k}, 1)
        error('RDCL_MorphImages all grid sizes must be equal.');
    end

    if size(varargin{k}, 3) ~= 1
        error('RDCL_MorphImages all images must be grayscale.');
    end
end

% Store images & scale
images(:, :, 1) = image1*morphProp(1);
images(:, :, 2) = image2*morphProp(2);

% Get extra images
for k = 1:length(varargin)
    images(:, :, k + 2) = varargin{k}*morphProp(k);
end

% Make them all doubles
images = double(images);

% How many images are there?
nImages = length(varargin) + 2;

% Morph them
morphedImage = round(sum(images, 3));
