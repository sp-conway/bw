function morphedDotPattern = RDCL_MorphDotPatterns(morphProp, dotPattern1, dotPattern2, varargin)

%
% Function: RDCL_MorphDotPatterns
%
% Morphs 2 or more dot patterns. All dot patterns must have the same number
% of dots and grid size. Simple linear combination.
%
% Arguments:
%   morphProp = 1xn vector of proportion of each dot pattern present in
%     morph. If doesn't sum to 1, will be rescaled. n = number of dot
%     patterns.
%   dotPattern1 = the first dot pattern.
%   dotPattern2 = the second dot pattern.
%   varargin:
%     dotPatternN = more dot patterns.
%
% Return:
%   morphedDotPattern = the morphed dot pattern.
%
% Example:
%   dotPattern1 = RDCL_GenerateDotPattern([500 500], 9);
%   dotPattern2 = RDCL_GenerateDotPattern([500 500], 9);
%   newDotPattern = RDCL_MorphDotPatterns([.5 .5], dotPattern1, dotPattern2);
%
% Notes:
%
% Andrew Cohen
% 6/25/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Make sure proportions sum to 1
morphProp = abs(morphProp);
morphProp = morphProp./sum(morphProp);

% Store dot patterns
dotPatterns{1} = dotPattern1;
dotPatterns{2} = dotPattern2;

% Get extra dot patterns
for k = 1:length(varargin)
    dotPatterns{2+k} = varargin{k};
end

% How many patterns are there?
nPatterns = length(varargin) + 2;

% All grid sizes must be the same
gridSize_x = dotPatterns{1}{2}(1);
gridSize_y = dotPatterns{1}{2}(2);

for i = 1:nPatterns
    if gridSize_x ~= dotPatterns{i}{2}(1)
        error('RDCL_MorphDotPatterns all grid sizes must be equal.');
    end

    if gridSize_y ~= dotPatterns{i}{2}(2)
        error('RDCL_MorphDotPatterns all grid sizes must be equal.');
    end
end


% All patterns must have the same number of dots
nDots = size(dotPatterns{1}{1}, 1);

for i = 1:nPatterns
    if nDots ~= size(dotPatterns{i}{1}, 1)
        error('RDCL_MorphDotPatterns all dot patterns must have the same number of dots.');
    end
end

% Morph them
morphedDotPattern{1} = repmat(0, [nDots, 2]);
for i = 1:nPatterns
    morphedDotPattern{1} = morphedDotPattern{1} + ...
        morphProp(i)*dotPatterns{i}{1};
end

morphedDotPattern{1} = round(morphedDotPattern{1});

% Copy over grid size
morphedDotPattern{2} = dotPattern1{2};
