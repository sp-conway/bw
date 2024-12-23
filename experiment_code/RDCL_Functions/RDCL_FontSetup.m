function [old_font_name old_font_number old_font_size old_font_style old_font_mode] = RDCL_FontSetup(varargin);

% function RDCL_FontSetup
%
% Sets up font defaults.
%
% Arguments:
%   screen_ptr = psyschtoolbox screen pointer.
%   varargin:
%     'TextFont' = font name
%     'TextSize' = font size
%     'TextStyle' = font style
%     'TextMode' = font mode (Screen('TextModes'))
%
% Returns:
%    old_font_name = font options before changes
%    old_font_number
%    old_font_size
%    old_font_style
%    old_font_mode
%
% Example:
%   [old_font_name old_font_number old_font_size old_font_style
%   old_font_mode] = RDCL_FontSetup(screen_ptr, 'Font', 'Arial',
%   'FontSize', 32, 'FontStyle', 0, 'TextMode', 'TextFill');
%
% Andrew Cohen
% 5/5/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

global screen_ptr

% Get old values
[old_font_name old_font_number] = Screen('TextFont', screen_ptr);
old_font_size = Screen('TextSize', screen_ptr);
old_font_style = Screen('TextStyle', screen_ptr);
old_font_mode = Screen('TextMode', screen_ptr);

k = 1;
while k < length(varargin)

    if strcmp('TextFont', varargin{k})
        font = varargin{k+1};
        Screen('TextFont', screen_ptr, font);
		k = k + 1;
    elseif strcmp('TextSize', varargin{k})
        font_size = varargin{k+1}(1);
        Screen('TextSize', screen_ptr, font_size);
		k = k + 1;
    elseif strcmp('TextStyle', varargin{k})
        font_style = varargin{k+1}(1);
        Screen('TextStyle', screen_ptr, font_style);
		k = k + 1;
    elseif strcmp('TextMode', varargin{k})
        font_mode = varargin{k+1}(1);
        Screen('TextMode', screen_ptr, font_mode);
		k = k + 1;
    else
        error('PACLAb_FontSetup: Incorrect option');
    end

	k = k + 1;

end
