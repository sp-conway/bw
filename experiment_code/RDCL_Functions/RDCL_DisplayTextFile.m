function RDCL_DisplayTextFile(text_file, varargin)

% function RDCL_DisplayTextFile
%
% Displays a text file.
%
% Arguments:
%    text_file = the name of the text file (including path).
%    varargin:
%      'Font' = the font name.
%      'FontSize' = the font size in pixels.
%      'bgColor' = the background color in RGB.
%      'Color' = the text color in RGB.
%      'LineHeight' = the height of a line in pixels.
%      'LeftMargin' = start of line from left of screen in pixels.
%      'TopMargin' = start of first line from top of screen in pixels.
%      'PauseTime' = how long to wait before accept button press in secs.
%      'InputDevice' = the number of the device that will get button press.
%      'LegalKey' = only a single key will exit the screen.
%      'StringSubs' = strings to substitute for '%s' in text file
%                     as in sprintf
%
% Returns:
%    none
%
% Example:
%   RDCL_DisplayTextFile('Text.txt', screen_ptr, 'Font', 'Arial', ...
%                          'FontSize', 24);
%
% Andrew Cohen
% 5/5/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
% 12/6/06 - Michael Ross - Modified to use new RDCL_GetResponse function
% 7/27/11 - alc - Moved constants from RDCL_Constants to this file.

global screen_ptr

% Defaults
text_font = 'Arial';
text_font_size = 32;
bg_color = RDCL_GetGray();
text_color = BlackIndex(screen_ptr);
text_line_height = -1;
text_line_left_margin = -1; % -1 = font size
text_line_top_margin = -1; % -1 = font size
pause_time = 3;
legal_input = [];
string_subs = [];
input_device = -1;

% if def_text_line_height = -1, line height set to def_text_line_percent*font size
def_text_line_percent = .2;

% Get user options
if nargin > 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'Font')
            text_font = varargin{k+1};
        elseif strcmp(option_text,'FontSize')
            text_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'bgColor')
            bg_color = varargin{k+1};
        elseif strcmp(option_text, 'Color')
            text_color = varargin{k+1};
        elseif strcmp(option_text, 'LineHeight')
            text_line_height = varargin{k+1}(1);
        elseif strcmp(option_text, 'LeftMargin')
            text_line_left_margin = varargin{k+1}(1);
        elseif strcmp(option_text, 'TopMargin')
            text_line_top_margin = varargin{k+1}(1);
        elseif strcmp(option_text, 'PauseTime')
            pause_time = varargin{k+1}(1);
        elseif strcmp(option_text, 'InputDevice')
            input_device = varargin{k+1}(1);
        elseif strcmp(option_text, 'LegalKey')
            legal_input = varargin{k+1}(1);
        elseif strcmp(option_text, 'StringSubs')
            string_subs = varargin{k+1};
        else
            error('RDCL_DisplayTextFile: Incorrect option');
        end

    end
end

% More default handling
if text_line_height == -1
    text_line_height = text_font_size + text_font_size*def_text_line_percent;
end
if text_line_left_margin == -1
    text_line_left_margin = text_font_size;
end
if text_line_top_margin == -1
    text_line_top_margin = text_font_size;
end

% Set options
[old_font_name old_font_number] = Screen('TextFont', screen_ptr, text_font);
old_text_size = Screen('TextSize', screen_ptr, text_font_size);

% Show the file
Screen('FillRect', screen_ptr, bg_color);
text_file_ptr = fopen(text_file);
line_num = 0;
while 1
    tmp_text = fgetl(text_file_ptr);
    if ~ischar(tmp_text)
        break
    else
        if ~isempty(string_subs) && ~isempty(strfind(tmp_text, '%s'))
            tmp_text = sprintf(tmp_text, string_subs{1});
            string_subs = string_subs(2:end);
        end

        Screen('DrawText', screen_ptr, tmp_text, text_line_left_margin, ...
               text_line_top_margin + line_num*text_line_height, text_color);
    end
    line_num = line_num + 1;
end
fclose(text_file_ptr);
Screen('Flip', screen_ptr);
pause(pause_time);
if isempty(legal_input)
	if input_device ~= -1
    	RDCL_WaitForKeyPress(input_device);
	else
	  	RDCL_WaitForKeyPress();
	end
else
	if input_device ~= -1
    	RDCL_GetResponse({legal_input}, input_device);
	else
		RDCL_GetResponse({legal_input});
	end
end

% Reset font
Screen('TextFont', screen_ptr, old_font_name);
Screen('TextSize', screen_ptr, old_text_size);
