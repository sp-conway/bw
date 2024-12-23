
function RDCL_ScreenCapture(file_name, varargin)

% function RDCL_ScreenCapture
%
% Captures and saves a copy of the current screen.
%
% Arguments:
%    file_name = the name of the saved file (including path). The name
%    should not include the file type extension.
%    varargin:
%      'FileType' = the type of file to save, e.g., 'bmp' or 'tiff'.
%        bmp is default. Not all types work.
%
% Returns:
%    none
%
% Example:
%   RDCL_ScreenCapture('tmp');
%   RDCL_ScreenCapture('tmp', 'FileType', 'tiff');
%
% Andrew Cohen
% 2/11/09
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none

global screen_ptr

% Constants
RDCL_Constants;

% Defaults
file_type = def_screen_capture_type;

% Get user options
if nargin > 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'FileType')
            file_type = varargin{k+1};
        else
            error('Incorrect DisplayTextFile option');
        end

    end
end

imageArray = Screen('GetImage', screen_ptr);
imwrite(imageArray, [file_name '.' file_type], file_type)
