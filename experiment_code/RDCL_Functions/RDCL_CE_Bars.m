function [] = RDCL_CE_Bars(bar_sizes, varargin)

% Displays stimuli for context effects as bars
%
%  bar_sizes = a matrix the bar sizes. rows = alts. cols = dims. Must be specified
%
%  varargin:
%
%    bar_loc: [X Y] in pixels or 'center'. Loction of center of display. 'center' = default.
%    bar_dir: 'v' or 'h'. Vertical or horizontal bars. 'h' = default. ONLY 'h' CURRENTLY IMPLEMENTED.
%    bar_width: 0-... in pixels. The width of the bars. 50 = default.
%    bar_dist: 0-... in pixels. Distance between bars within a group. 10 = default.
%    bar_group_dist: 0-... in pixels. Distance between each group of bars. 50 = default.
%    bar_col: [r g b] bar color. All bars are the same color. [0 0 255] = default.
%    bar_start_line: boolean. true = default.
%
%    tick_pos: vector of pixels. Where to put the axis ticks. NA = default. Largest tick_pos defines max bar size. If NA, uses max bar_sizes.
%    tick_len: int. length of a tick. 5 = default.
%    tick_lab: cell array of strings. Labels for each tick. {} = default.
%    tick_lab_col: [r g b] axis color. Also the color of the axis line. [255 255 255] = default.
%    tick_lab_font_size: int. 14 = default.
%    tick_title: string. tick title label. '' = default.
%    tick_title_col: [r g b] tick title color [255 255 255] = default.
%    tick_title_font_size: int. 18 = default.
%    tick_hide: boolean. false = default.
%
%    alt_lab: cell array of strings. Labels for each alternative. {A1 A2 A3...} = default.
%    alt_lab_col: [r g b] label color for alternatives. [255 255 255] = default.
%    alt_lab_font_size: int. 18 = default.
%    alt_title: string. Alternative label. NA = default.
%    alt_title_col: [r g b] label color for alternatives. [255 255 255] = default.
%    alt_title_font_size: int. 24 = default.
%
%    dim_lab: cell array of strings. Labels for each dimension. {D1 D2 D3...} = default.
%    dim_lab_col: [r g b] label color for dimensions. [255 255 255] = default
%    dim_lab_font_size: int. 18 = default.
%    dim_title: string. Dimension label. NA = default.
%    dim_title_col: [r g b] label color for dimensions. [255 255 255] = default
%    dim_title_font_size: int. 24 = default.
%
% Returns
%   Nothing
%
% Notes:
%
% Example:
%
%
% Andrew Cohen
% 10/4/19
% Copyright (c) 2019 Andrew L. Cohen
%
% Change notes
%

% Constants -------------------------------------------

global screen_ptr;

% Some counts
n_alt = size(bar_sizes)(1);
n_dim = size(bar_sizes)(2);

% The screen size and center
screen_rect = Screen('Rect', screen_ptr);
screen_cent = [round(screen_rect(3)/2) round(screen_rect(4)/2)];

% Defaults --------------------

bar_loc = default_bar_loc = screen_cent;
bar_dir = default_bar_dir =  'h';
bar_width = default_bar_width =  50;
bar_dist = default_bar_dist =  10;
bar_group_dist = default_bar_group_dist =  50;
bar_col = default_bar_col =  [0 0 255];
bar_start_line = default_bar_start_line = true; 

tick_pos = default_tick_pos =  NA;
tick_len = default_tick_len = 5;
tick_lab = default_tick_lab =  {};
tick_lab_col = default_tick_lab_col =  [255 255 255];
tick_lab_font_size = default_tick_lab_font_size = 14;
tick_title = default_tick_title = '';
tick_title_col = default_tick_title_col = [255 255 255];
tick_title_font_size = default_tick_title_font_size = 18;
tick_hide = default_tick_hide = false;

alt_lab = default_alt_lab = strsplit(strtrim(sprintf('A%d ', 1:n_alt)), del=' ');
alt_lab_col = default_alt_lab_col = [255 255 255];
alt_lab_font_size = default_alt_lab_font_size = 18;
alt_title = default_alt_title = '';
alt_title_col = default_alt_title_col = [255 255 255];
alt_title_font_size = default_alt_title_font_size = 24;

dim_lab = default_dim_lab = strsplit(strtrim(sprintf('D%d ', 1:n_dim)), del=' ');
dim_lab_col = default_dim_lab_col = [255 255 255];
dim_lab_font_size = default_dim_lab_font_size = 18;
dim_title = default_dim_title = '';
dim_title_col = default_dim_title_col = [255 255 255];
dim_title_font_size = default_dim_title_font_size = 24;

% Get user options -----------------
if nargin > 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'bar_loc')
          bar_loc = varargin{k+1};
        elseif strcmp(option_text, 'bar_dir')
          bar_dir = varargin{k+1};
        elseif strcmp(option_text, 'bar_width')
          bar_width = varargin{k+1}(1);
        elseif strcmp(option_text, 'bar_dist')
          bar_dist = varargin{k+1}(1);
        elseif strcmp(option_text, 'bar_group_dist')
          bar_group_dist = varargin{k+1}(1);
        elseif strcmp(option_text, 'bar_col')
          bar_col = varargin{k+1};
        elseif strcmp(option_text, 'tick_pos')
          tick_pos = varargin{k+1};
        elseif strcmp(option_text, 'tick_len')
          tick_len = varargin{k+1}(1);
        elseif strcmp(option_text, 'tick_lab')
          tick_lab = varargin{k+1};
        elseif strcmp(option_text, 'tick_lab_col')
          tick_lab_col = varargin{k+1};
        elseif strcmp(option_text, 'tick_lab_font_size')
          tick_lab_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'tick_title')
          tick_title = varargin{k+1};
        elseif strcmp(option_text, 'tick_title_col')
          tick_title_col = varargin{k+1};
        elseif strcmp(option_text, 'tick_title_font_size')
          tick_title_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'tick_hide')
          tick_title_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'alt_lab')
          alt_lab = varargin{k+1};
        elseif strcmp(option_text, 'alt_lab_col')
          alt_lab_col = varargin{k+1};
        elseif strcmp(option_text, 'alt_lab_font_size')
          alt_lab_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'alt_title')
          alt_title = varargin{k+1};
        elseif strcmp(option_text, 'alt_title_col')
          alt_title_col = varargin{k+1};
        elseif strcmp(option_text, 'alt_title_font_size')
          alt_title_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'dim_lab')
          dim_lab = varargin{k+1};
        elseif strcmp(option_text, 'dim_lab_col')
          dim_lab_col = varargin{k+1};
        elseif strcmp(option_text, 'dim_lab_font_size')
          dim_lab_font_size = varargin{k+1}(1);
        elseif strcmp(option_text, 'dim_title')
          dim_title = varargin{k+1};
        elseif strcmp(option_text, 'dim_title_col')
          dim_title_col = varargin{k+1};
        elseif strcmp(option_text, 'dim_title_font_size')
          dim_title_font_size = varargin{k+1}(1);
        else
            error('RDCL_DisplayTextFile: Incorrect option');
        end

    end
end

% Get some sizes -------------------

% Max bar size
if (isna(tick_pos))
  max_bar_size = max(bar_sizes(:));
else
  max_bar_size = max(tick_pos);
end

% Max tick lab size
max_tick_lab_size = 0;
for (i = 1:length(tick_lab))

  woff = Screen('OpenOffscreenWindow', screen_ptr, [], [0 0 2*tick_lab_font_size*length(tick_lab{i}) 2*tick_lab_font_size]);
  Screen(woff, 'TextSize', tick_lab_font_size);
  bounds=TextBounds(woff, tick_lab{i}, 1);
  if (bounds(3)-bounds(1) > max_tick_lab_size)
    max_tick_lab_size = bounds(3)-bounds(1);
  end

end

% Max alt lab size
max_alt_lab_size = 0;
for (i = 1:n_alt)
  
  woff = Screen('OpenOffscreenWindow', screen_ptr, [], [0 0 2*alt_lab_font_size*length(alt_lab{i}) 2*alt_lab_font_size]);
  Screen(woff, 'TextSize', alt_lab_font_size);
  bounds=TextBounds(woff, alt_lab{i}, 1);
  if (bounds(3)-bounds(1) > max_alt_lab_size)
    max_alt_lab_size = bounds(3)-bounds(1);
  end
  
end

% Alt title size (rotated)
if (strcmp(alt_title, ''))
  alt_title_size = 0;
else
  woff = Screen('OpenOffscreenWindow', screen_ptr, [], [0 0 2*alt_title_font_size*length(alt_title) 2*alt_title_font_size]);
  Screen(woff, 'TextSize', alt_title_font_size);
  bounds=TextBounds(woff, alt_title, 1);
  alt_title_size = bounds(4) - bounds(2);
end

% Figure out some positions ------------------

% bar, title, label, tick positions
if (bar_dir == 'h')

  % min x start for bars
  bar_x_min = bar_loc(1) - round((max_bar_size*n_dim + (n_dim-1)*bar_group_dist)/2);

  % min y start for bars
  bar_y_min = bar_loc(2) - round((bar_width*n_alt + (n_alt-1)*bar_dist)/2);

  % x-position for bar starts
  bar_x = repmat(bar_x_min, [1 n_dim]) + (0:(max_bar_size+bar_group_dist):((n_dim-1)*(max_bar_size+bar_group_dist)));

  % y-position for bar starts
  bar_y = repmat(bar_y_min, [1 n_alt]) + (0:(bar_width+bar_dist):((n_alt-1)*(bar_width+bar_dist)));
        RDCL_DrawText(1, 1, ' ', 'Col', tick_lab_col); % FONT SIZE?

  % x-position for alt labels
  alt_lab_x = bar_x_min - bar_group_dist - max_alt_lab_size;

  % y-position for alt labels
  alt_lab_y = bar_y + round(bar_width/2);

  % x-position for alt title
  alt_title_x = alt_lab_x - 2*alt_title_size;

  % y-position for alt title
  alt_title_y = bar_loc(2);

  % x-position for dim labels
  dim_lab_x = bar_x + round(max_bar_size)/2;

  % y-position for dim labels
  dim_lab_y = bar_y_min - bar_dist;

  % x-position for dim title
  dim_title_x = bar_loc(1);

  % y-position for dim title
  dim_title_y = dim_lab_y - (bar_width/2 + bar_dist);

  % x-position for axis
  tick_axis_x = bar_x;

  % y-position for axis
  tick_axis_y = bar_y(n_alt) + bar_width + bar_dist;

  % x-position for ticks
  tick_x = tick_pos;

  % y-position for ticks
  tick_y = tick_axis_y;
  
  % x-position for tick labels
  tick_lab_x = tick_x;

  % y-position for tick labels
  tick_lab_y = tick_axis_y + 1.5*tick_len + 1.5*max_tick_lab_size;

  % x-position for tick title
  tick_title_x = bar_x + round(max_bar_size/2);

  % y-position for tick title
  tick_title_y = tick_lab_y + bar_dist;

elseif (bar_dir == 'v')

  % NOT IMPLEMENTED

end

% Kludge -----------------
% Without this, DrawFormattedText2 wouldn't work
RDCL_DrawText(1, 1, ' ', 'Col', tick_lab_col); % FONT SIZE?
% B/c of use of offscreen buffer above? TextBounds?

% Axis --------------------------

tick_axis_col = tick_lab_col;
tick_col = tick_lab_col;

axis_pen_width = 1;

% Axis line
if (!isna(tick_pos))
  for (i = 1:n_dim)
    Screen('DrawLine', screen_ptr, tick_axis_col, ...
        tick_axis_x(i), tick_axis_y, tick_axis_x(i) + max_bar_size, tick_axis_y, axis_pen_width);
  end
end

n_ticks = length(tick_pos);

% Ticks and labels
if (!isna(tick_pos))
  for (i = 1:n_dim)
    for (j = 1:n_ticks)
      
      the_x = tick_x(j) + bar_x(i);
      start_y = tick_y - tick_len;
      end_y = tick_y + tick_len;
      Screen('DrawLine', screen_ptr, tick_col, ...
          the_x, start_y, the_x, end_y, axis_pen_width);
          
      % Only aligns well if you don't use descenders
      if (!isempty(tick_lab))
        tick_lab_tmp = ['<size=', num2str(tick_lab_font_size), '>', tick_lab{j}];
        DrawFormattedText2(tick_lab_tmp, 
          'win', screen_ptr, 
          'sx', the_x, 'sy', tick_lab_y,
          'xalign', 'center', 'yalign', 'bottom',
          'baseColor', tick_lab_col);
      end
          
    end
  end
end

% Axis title
if (!strcmp(tick_title, ''))
  for (i = 1:n_dim)
    tick_title_tmp = ['<size=', num2str(tick_title_font_size), '>', tick_title];
    DrawFormattedText2(tick_title_tmp, 
      'win', screen_ptr, 
      'sx', tick_title_x(i), 'sy', tick_title_y,
      'xalign', 'center', 'yalign', 'top',
      'baseColor', tick_lab_col);
    %RDCL_DrawText(tick_title_x(i), tick_title_y, tick_title, 'Col', dim_lab_col); % FONT SIZE?
  end
end

% Labels --------------------------

% Alternative title
if (!strcmp(alt_title, ''))
  DrawFormattedText2(alt_title, 
    'win', screen_ptr, 
    'sx', alt_title_x, 'sy', alt_title_y,
    'xalign', 'center', 'yalign', 'center',
    'baseColor', alt_title_col,
    'transform', {'rotate', -90});
end

% Alternative labels
if (!isempty(alt_lab))
  for (i = 1:n_alt)
    alt_lab_tmp = ['<size=', num2str(alt_lab_font_size), '>', alt_lab{i}];
    DrawFormattedText2(alt_lab_tmp, 
      'win', screen_ptr, 
      'sx', alt_lab_x, 'sy', alt_lab_y(i),
      'xalign', 'left', 'yalign', 'center',
      'baseColor', alt_lab_col);
  end
end

% Dimension title
if (!strcmp(dim_title, ''))
  dim_title_tmp = ['<size=', num2str(dim_title_font_size), '>', dim_title];
  DrawFormattedText2(dim_title_tmp, 
    'win', screen_ptr, 
    'sx', dim_title_x, 'sy', dim_title_y,
    'xalign', 'center', 'yalign', 'bottom',
    'baseColor', dim_title_col);
end

% Dimension labels
if (!isempty(dim_lab))
  for (i = 1:n_dim)
      dim_lab_tmp = ['<size=', num2str(dim_lab_font_size), '>', dim_lab{i}];
      DrawFormattedText2(dim_lab_tmp, 
        'win', screen_ptr, 
        'sx', dim_lab_x(i), 'sy', dim_lab_y,
        'xalign', 'center', 'yalign', 'bottom',
        'baseColor', dim_lab_col);
  end
end

% Bars --------------------------

% Bars
for (i = 1:n_alt)
  for (j = 1:n_dim)
    
    bar_x_start = bar_x(j);
    bar_y_start = bar_y(i);
    bar_x_end = bar_x_start + bar_sizes(i, j);
    bar_y_end = bar_y_start + bar_width;
    bar = [bar_x_start bar_y_start bar_x_end bar_y_end];
    Screen('FillRect', screen_ptr, bar_col, bar);
    
  end
end

% Start line
bar_start_line_col = tick_col;

for (i = 1:n_dim)
  Screen('DrawLine', screen_ptr, bar_start_line_col, ...
      bar_x(i), bar_y(1), bar_x(i), bar_y(n_alt)+bar_width, axis_pen_width);
end

endfunction
