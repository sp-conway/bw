function [button_pressed button_accepted] = RDCL_Button(button_pressed, varargin)

% function RDCL_Button
%
% Implements a slider bar.
%
% Arguments:
%    varargin:
%      'ButtonColor' = Button color.
%      'ButtonColorPressed' = Button color when pressed.
%      'FrameColor' = Button frame color.
%      'FrameColorPressed' = Button frame color when pressed.
%      'FrameLineWidth' = Width of frame.
%      'ButtonText' = Text in the button.
%      'ButtonTextPressed' = Text in the button when pressed.
%      'ButtonTextColor' = The color of the text.
%      'ButtonTextColorPressed' = The color of the text when pressed.
%      'Font' = Font.
%      'FontSize' = Font size.
%      'ButtonHeight' = Button height.
%      'ButtonWidth' = ButtonWidth.
%      'CenterX' = Button center x location.
%      'CenterY' = Button center y location.
%      'ButtonImage' = Image to show instead of a rectangle.
%      'ButtonImagePressed' = Image to show instead of a rectangle when pressed.
%      'ButtonImageFile' = Image file to show instead of a rectangle.
%      'ButtonImageFilePressed' = Image file to show instead of a rectangle when pressed.
%      'AcceptOnUp' = Accept on button up (true) or down (false)
%
% Returns:
%    button_pressed = true if pressed, otherwise false.
%    button_accepted = true if the button was accepted, false otherwise.
%      AcceptOnUp option determines behavior.
%
% Example:
%
%    done = false;
%    button_pressed = false;
%    while (~done)
%      [button_pressed button_accepted] = RDCL_Button(button_pressed);
%      Screen('Flip', screen_ptr);
%      if button_accepted
%        done = true;
%      endif
%    end
%
%   Could also use:
%
%     button_image = imread('../RDCL_Resources/car_base.png');
%     [button_pressed button_accepted] = RDCL_Button(button_pressed, 'ButtonImage', button_image);
%
%     [button_pressed button_accepted] = RDCL_Button(button_pressed, ...
%       'ButtonImageFile', '../RDCL_Resources/car_base.png',
%       'ButtonImageFilePressed', '../RDCL_Resources/rdcl_logo.png', ...
%       'AcceptOnUp', true,
%       'FrameColorPressed', [244 0 0], ...
%       'FrameLineWidth', 5);
%
% Andrew Cohen
% 5/25/19
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%   For RT purposes, accepts click when button down by default, not when button up.
%   If an image is used, it is scaled to the button width and height.
%
% Known Bugs:
%   none
%
% Change History:
%

global screen_ptr;

% Get info
[screen_width, screen_height]=Screen('WindowSize', screen_ptr);

prev_button_pressed = button_pressed;

% Defaults
button_color = RDCL_GetGray();
frame_color = [255 255 255];
button_color_pressed = [255 255 255];
frame_color_pressed = [255 255 255];
frame_line_width = 1;

button_text = 'Yes';
button_text_pressed = [];
button_text_color = [0 0 0];
button_text_color_pressed = [];

text_font = 'Arial';
text_font_size = 24;

button_height = 50;
button_width = 100;

center_x = round(screen_width/2);
center_y = round(screen_height/2);

button_image = [];
button_image_pressed = [];
button_image_file = [];
button_image_file_pressed = [];

accept_on_up = false;

% Get user options
if nargin > 1
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'ButtonColor')
            button_color = varargin{k+1};
        elseif strcmp(option_text, 'FrameColor')
            frame_color = varargin{k+1};
        elseif strcmp(option_text, 'FrameLineWidth')
            frame_line_width = varargin{k+1};
        elseif strcmp(option_text, 'ButtonColorPressed')
            button_color_pressed = varargin{k+1};
        elseif strcmp(option_text, 'FrameColorPressed')
            frame_color_pressed = varargin{k+1};
        elseif strcmp(option_text, 'ButtonText')
            button_text = varargin{k+1};
        elseif strcmp(option_text, 'ButtonTextPressed')
            button_text_pressed = varargin{k+1};
        elseif strcmp(option_text, 'ButtonTextColor')
            button_text_color = varargin{k+1};
        elseif strcmp(option_text, 'ButtonTextColorPressed')
            button_text_color_pressed = varargin{k+1};
        elseif strcmp(option_text, 'Font')
            text_font = varargin{k+1};
        elseif strcmp(option_text, 'FontSize')
            text_font_size = varargin{k+1};
        elseif strcmp(option_text, 'ButtonHeight')
            button_height = varargin{k+1};
        elseif strcmp(option_text, 'ButtonWidth')
            button_width = varargin{k+1};
        elseif strcmp(option_text, 'CenterX')
            center_x = varargin{k+1};
        elseif strcmp(option_text, 'CenterY')
            center_y = varargin{k+1};
        elseif strcmp(option_text, 'ButtonImage')
            button_image = varargin{k+1};
        elseif strcmp(option_text, 'ButtonImagePressed')
            button_image_pressed = varargin{k+1};
        elseif strcmp(option_text, 'ButtonImageFile')
            button_image_file = varargin{k+1};
        elseif strcmp(option_text, 'ButtonImageFilePressed')
            button_image_file_pressed = varargin{k+1};
        elseif strcmp(option_text, 'AcceptOnUp')
            accept_on_up = varargin{k+1};
        else
            error('RDCL_Button: Incorrect option');
        end

    end
end

% Set font
[old_font_name old_font_number] = Screen('TextFont', screen_ptr, text_font);
old_text_size = Screen('TextSize', screen_ptr, text_font_size);

% Some other options
if isempty(button_text_pressed)
  button_text_pressed = button_text;
endif

if isempty(button_text_color_pressed)
  button_text_color_pressed = button_text_color;
endif

if isempty(button_image_pressed)
  button_image_pressed = button_image;
endif

if isempty(button_image_file_pressed)
  button_image_file_pressed = button_image_file;
endif

% Button rectangle
button_rect = [round(center_x - button_width/2) ...
               round(center_y - button_height/2) ...
               round(center_x + button_width/2) ...
               round(center_y + button_height/2)];
frame_rect = button_rect + ...
  [-frame_line_width -frame_line_width frame_line_width frame_line_width];

% Mouse input
[mouse_x mouse_y mouse_buttons] = GetMouse();

if (mouse_x < button_rect(3) && ...
    mouse_x > button_rect(1) && ...
    mouse_y < button_rect(4) && ...
    mouse_y > button_rect(2))

    if mouse_buttons(1) == 1

      button_pressed = true;

      if accept_on_up == false
        button_accepted = true;
      else
        button_accepted = false;
      endif

    else

      button_pressed = false;

      if accept_on_up == true && prev_button_pressed == true
        button_accepted = true;
      else
        button_accepted = false;
      endif

    endif

else

  button_pressed = false;
  button_accepted = false;

endif

% Draw it
if !isempty(button_image)

  if (button_pressed == true)
    Screen('PutImage', screen_ptr, button_image_pressed, button_rect);
    Screen('FrameRect', screen_ptr, frame_color_pressed, frame_rect, frame_line_width);
  else
    Screen('PutImage', screen_ptr, button_image, button_rect);
    Screen('FrameRect', screen_ptr, frame_color, frame_rect, frame_line_width);
  end

elseif !isempty(button_image_file)

  if (button_pressed == true)
    button_image = imread(button_image_file_pressed);
    Screen('PutImage', screen_ptr, button_image, button_rect);
    Screen('FrameRect', screen_ptr, frame_color_pressed, frame_rect, frame_line_width);
  else
    button_image = imread(button_image_file);
    Screen('PutImage', screen_ptr, button_image, button_rect);
    Screen('FrameRect', screen_ptr, frame_color, frame_rect, frame_line_width);
  endif

else

  if (button_pressed == true)

    Screen('FillRect', screen_ptr, button_color_pressed, button_rect);
    Screen('FrameRect', screen_ptr, frame_color_pressed, frame_rect, frame_line_width);
    [normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = ...
      Screen('TextBounds', screen_ptr, button_text_pressed);
    text_x = round(button_rect(1) + (button_rect(3) - button_rect(1))/2 - normBoundsRect(3)/2);
    text_y = round(button_rect(2) + (button_rect(4) - button_rect(2))/2 - normBoundsRect(4)/2);
    Screen('DrawText', screen_ptr, button_text_pressed, text_x, text_y, button_text_color_pressed);

  else

    Screen('FillRect', screen_ptr, button_color, button_rect);
    Screen('FrameRect', screen_ptr, frame_color, frame_rect, frame_line_width);
    [normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = ...
      Screen('TextBounds', screen_ptr, button_text);
    text_x = round(button_rect(1) + (button_rect(3) - button_rect(1))/2 - normBoundsRect(3)/2);
    text_y = round(button_rect(2) + (button_rect(4) - button_rect(2))/2 - normBoundsRect(4)/2);
    Screen('DrawText', screen_ptr, button_text, text_x, text_y, button_text_color);

  endif

endif

% Reset font
Screen('TextFont', screen_ptr, old_font_name);
Screen('TextSize', screen_ptr, old_text_size);
