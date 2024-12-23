function bar_val = RDCL_SliderBar(bar_val, varargin)

% function RDCL_SliderBar
%
% Implements a slider bar.
%
% Arguments:
%    bar_val = the start value of the bar.
%    varargin:
%      'BarMin' = Minimum bar value.
%      'BarMax' = Maximum var value.
%      'BarWidth' = Bar width in pixels.
%      'BarHeight' = Bar height in pixels.
%      'CenterX' = Bar center x location.
%      'CenterY' = Bar center y location.
%      'BarColor' = Bar color. [R B G] 0-255.
%      'BarBGColor' = Bar background color. [R B G] 0-255.
%      'FrameColor' = Frame color.
%      'TickVals' = Vector of tick locations.
%      'TickLen' = Length of the tick
%      'TickText' = Vector of tick strings, must match length of TickVals.
%         If none, value of TickVals. Format {'a', 'b', 'c'}.
%      'MouseBufferX' = Pixels at min and max vals that still record mouse location.
%      'MouseBufferY' = Pixels at min and max vals that still record mouse location.
%      'TextBuffer' = Separation of text and frame.
%      'BarBuffer' = Separation of bar and frame.
%      'LineWidth' = Width of frame.
%      'ShowCurrVal' = Show current value true or false.
%      'ShowTicks' = Show ticks true or false.
%      'CurrValPrecision' = Current value precision in decimals (integer).
%      'Font' = Font.
%      'FontSize' = Font size.
%
% Returns:
%    bar_val = the current value of the bar.
%
% Examples:
%
%   done = false;
%   bar_val = 25;
%   while (~done)
%     bar_val = RDCL_SliderBar(bar_val);
%     Screen('Flip', screen_ptr);
%     % if ... done = true endif
%   end
%
%   done = false;
%   bar_vals = [25 50 25];
%   prev_bar_vals = bar_vals;
%   total = sum(bar_vals);
%   while (~done)
%
%     bar_vals(1) = RDCL_SliderBar(bar_vals(1), 'CenterY', 200);
%     bar_vals(2) = RDCL_SliderBar(bar_vals(2), 'CenterY', 400);
%     bar_vals(3) = RDCL_SliderBar(bar_vals(3), 'CenterY', 600);
%
%     changed = find (bar_vals != prev_bar_vals);
%     not_changed = find (bar_vals == prev_bar_vals);
%
%     if (!isempty(changed))
%       if (sum(bar_vals(not_changed))> 0)
%         bar_vals(not_changed) = ...
%           (total - bar_vals(changed))*bar_vals(not_changed)/sum(bar_vals(not_changed));
%       else
%         bar_vals(not_changed) = ...
%           (total - bar_vals(changed))/length(not_changed);
%       endif
%     endif
%
%     prev_bar_vals = bar_vals;
%
%     Screen('Flip', screen_ptr);
%
%   end
%
% Andrew Cohen
% 5/24/19
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%   Does not Flip the screen or collect the response.
%   This way other things can be drawn on the screen with the slider.
%   See the example for how to draw the slider. The example never breaks.
%     Put in some sort of response to stop it.
%
% Known Bugs:
%   none
%
% Change History:
%

global screen_ptr;

% Get info
[screen_width, screen_height]=Screen('WindowSize', screen_ptr);

% Defaults
bar_min_val = 0;
bar_max_val = 100;

bar_width_max = 200;
bar_height_max = 40;

bar_color_inner = RDCL_GetGray();
bar_bg_color_inner = [0 0 0];
bar_color_outer = RDCL_GetGray();

bar_center_x = round(screen_width/2);
bar_center_y = round(screen_height/2);

tick_locs = [bar_min_val RDCL_Round((bar_min_val + bar_max_val)/2, 2) bar_max_val];
tick_text = {};
tick_line_len = 5;

mouse_buffer_x = 20;
mouse_buffer_y = 20;
text_buffer = 5;
bar_buffer = 1;

line_width = 2;

show_curr_val = true;
show_ticks = true;
curr_val_precision = 0;

text_font = 'Arial';
text_font_size = 32;

% Get user options
if nargin > 1
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'BarMin')
            bar_min_val = varargin{k+1};
        elseif strcmp(option_text, 'BarMax')
            bar_max_val = varargin{k+1};
        elseif strcmp(option_text, 'BarWidth')
            bar_width_max = varargin{k+1};
        elseif strcmp(option_text, 'BarHeight')
            bar_height_max = varargin{k+1};
        elseif strcmp(option_text, 'CenterX')
            bar_center_x = varargin{k+1};
        elseif strcmp(option_text, 'CenterY')
            bar_center_y = varargin{k+1};
        elseif strcmp(option_text, 'BarColor')
            bar_color_inner = varargin{k+1};
        elseif strcmp(option_text, 'BarBGColor')
            bar_bg_color_inner = varargin{k+1};
        elseif strcmp(option_text, 'FrameColor')
            bar_color_outer = varargin{k+1};
        elseif strcmp(option_text, 'TickVals')
            tick_locs = varargin{k+1};
        elseif strcmp(option_text, 'TickLen')
            tick_line_len = varargin{k+1};
        elseif strcmp(option_text, 'TickText')
            tick_text = varargin{k+1};
        elseif strcmp(option_text, 'MouseBufferX')
            mouse_buffer_x = varargin{k+1};
        elseif strcmp(option_text, 'MouseBufferY')
            mouse_buffer_y = varargin{k+1};
        elseif strcmp(option_text, 'TextBuffer')
            text_buffer = varargin{k+1};
        elseif strcmp(option_text, 'BarBuffer')
            bar_buffer = varargin{k+1};
        elseif strcmp(option_text, 'LineWidth')
            line_width = varargin{k+1};
        elseif strcmp(option_text, 'ShowCurrVal')
            show_curr_val = varargin{k+1};
        elseif strcmp(option_text, 'ShowTicks')
            show_ticks = varargin{k+1};
        elseif strcmp(option_text, 'CurrValPrecision')
            curr_val_precision = varargin{k+1};
        elseif strcmp(option_text, 'Font')
            text_font = varargin{k+1};
        elseif strcmp(option_text, 'FontSize')
            text_font_size = varargin{k+1};
        else
            error('RDCL_SliderBar: Incorrect option');
        end

    end
end

[old_font_name old_font_number] = Screen('TextFont', screen_ptr, text_font);
old_text_size = Screen('TextSize', screen_ptr, text_font_size);

bar_height = bar_height_max;

bar_inner_min_x = round(bar_center_x - bar_width_max/2);
bar_inner_max_x = round(bar_center_x + bar_width_max/2);
bar_inner_min_y = round(bar_center_y - bar_height/2);
bar_inner_max_y = round(bar_center_y + bar_height/2);

if length(tick_text) == 0
  for i = 1:length(tick_locs)
    tick_text{i} = num2str(tick_locs(i));
  endfor
endif

% Get user mouse input
[mouse_x mouse_y mouse_buttons] = GetMouse();

if (mouse_x < bar_inner_max_x && ...
    mouse_x > bar_inner_min_x && ...
    mouse_y < bar_inner_max_y + mouse_buffer_y && ...
    mouse_y > bar_inner_min_y - mouse_buffer_y &&
    mouse_buttons(1) == 1)

    bar_val = bar_min_val + ...
      (bar_max_val - bar_min_val) ...
      *(mouse_x - bar_inner_min_x)/(bar_inner_max_x - bar_inner_min_x);

elseif (mouse_x < bar_inner_min_x && ...
        mouse_x > bar_inner_min_x - mouse_buffer_x && ...
        mouse_y < bar_inner_max_y + mouse_buffer_y  && ...
        mouse_y > bar_inner_min_y - mouse_buffer_y &&
        mouse_buttons(1) == 1)

    bar_val = bar_min_val;

elseif (mouse_x > bar_inner_max_x && ...
        mouse_x < bar_inner_max_x + mouse_buffer_x && ...
        mouse_y < bar_inner_max_y + mouse_buffer_y  && ...
        mouse_y > bar_inner_min_y - mouse_buffer_y &&
        mouse_buttons(1) == 1)

    bar_val = bar_max_val;

endif

% Make sure value is within bounds
bar_val = max(min(bar_val, bar_max_val), bar_min_val);

bar_width = round( bar_width_max*(bar_val - bar_min_val)/(bar_max_val - bar_min_val) );

% Outer rectangle
bar_rect_outer = [round(bar_center_x - bar_width_max/2 - line_width - bar_buffer) ...
                  round(bar_center_y - bar_height_max/2 - line_width - bar_buffer) ...
                  round(bar_center_x + bar_width_max/2 + line_width + bar_buffer) ...
                  round(bar_center_y + bar_height_max/2 + line_width + bar_buffer)];
Screen('FrameRect', screen_ptr, bar_color_outer, bar_rect_outer, line_width);

% Inner rectangle background
bar_rect_inner_bg = [round(bar_center_x - bar_width_max/2) ...
                     round(bar_center_y - bar_height_max/2) ...
                     round(bar_center_x + bar_width_max/2) ...
                     round(bar_center_y + bar_height_max/2)];
Screen('FillRect', screen_ptr, bar_bg_color_inner, bar_rect_inner_bg);

% Inner rectangle
bar_rect_inner = [round(bar_center_x - bar_width_max/2) ...
                  round(bar_center_y - bar_height/2) ...
                  round(bar_center_x - bar_width_max/2 + bar_width) ...
                  round(bar_center_y + bar_height/2)];
Screen('FillRect', screen_ptr, bar_color_inner, bar_rect_inner);

% Ticks
if (show_ticks)

  % Find max tick text height
  the_text = 'lj';
  [normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = ...
    Screen('TextBounds', screen_ptr, the_text);
    max_text_height = normBoundsRect(4);

  % Tick lines and labels
  for (i = 1:length(tick_locs))

    line_x = round( bar_inner_min_x + ...
      (tick_locs(i)-bar_min_val)/(bar_max_val-bar_min_val)*...
      (bar_inner_max_x - bar_inner_min_x) );
    Screen('DrawLine', screen_ptr, bar_color_outer, ...
      line_x, bar_rect_outer(4), ...
      line_x, bar_rect_outer(4) + tick_line_len, ...
      line_width)

    the_text = tick_text{i};
    [normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = ...
      Screen('TextBounds', screen_ptr, the_text);
    text_x = round( bar_rect_outer(1) + ...
      (tick_locs(i)-bar_min_val)/(bar_max_val-bar_min_val)*...
      (bar_rect_outer(3) - bar_rect_outer(1)) - ...
      normBoundsRect(3)/2);
    text_y = bar_rect_outer(4) + text_buffer + tick_line_len + max_text_height;
    Screen('DrawText', screen_ptr, the_text, text_x, text_y, bar_color_outer, [], 1);

  endfor

endif

% Current value
if show_curr_val

  line_x = bar_rect_inner(3);
  Screen('DrawLine', screen_ptr, bar_color_outer, ...
    line_x, bar_rect_outer(2), ...
    line_x, bar_rect_outer(2) - tick_line_len, ...
    line_width)

  the_text = sprintf(sprintf('%%.%df', curr_val_precision), bar_val);
  [normBoundsRect, offsetBoundsRect, textHeight, xAdvance] = ...
    Screen('TextBounds', screen_ptr, the_text);
  text_x = bar_rect_inner(3) - round(normBoundsRect(3)/2);
  text_y = round(bar_rect_outer(2) - normBoundsRect(4) - text_buffer - tick_line_len);
  Screen('DrawText', screen_ptr, the_text, text_x, text_y, bar_color_outer);

endif

% Reset font
Screen('TextFont', screen_ptr, old_font_name);
Screen('TextSize', screen_ptr, old_text_size);
