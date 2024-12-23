function spoke_locs = RDCL_SpokeDisplay(spoke_options)

% Displays a spoke graph
%
%  spoke_options = a struct with the arguments to make the bars
%     .dim_1_level: 0-1
%     .dim_2_level: 0-1
%     .dim_3_level: 0-1
%     .dim_4_level: 0-1
%     .axis_size: 0 = smallest, 1 = takes up whole screen
%         Uses the smaller of the x and y dimensions for max size
%     .axis_color: [r g b] color
%     .labels: cell of 4 labels (top, right, bottom, left)
%     .title: overall title
%
%  spoke_locs = a struct with the locations of the axes and lines
%    .dim_12_line (line from dim 1 to dim 2)
%    .dim_23_line
%    .dim_34_line
%    .dim_41_line
%    .vert_axis
%    .horiz_axis
%
% Example:
%
%   spoke_options.dim_1_level = .5;
%   spoke_options.dim_2_level = .5;
%   spoke_options.dim_3_level = .5;
%   spoke_options.dim_4_level = .5;
%   spoke_options.axis_size = .5;
%   spoke_options.axis_color = [255 255 255];
%   spoke_options.labels = {'10K', 'TOT. MILES', 'AGE', 'MOT.'};
%   spoke_options.title = 'Runner = 1';
%
%   Spoke_Display(spoke_options)
%
% Notes:
%   dim_1 is top, moving clockwise
%
% Andrew Cohen
% 1/14/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change log:
%   5/24/19 alc added to RDCL functions

% Constants -------------------------------------------

global screen_ptr;

% Defaults
default_font = 'Arial';
default_font_size = 26;

default_num_ticks = 10; % Not incl. origin
default_tick_size_ratio = .01;

pen_width = 1;

text_label_buffer = 1.5; % of vert extent

% Get options
dim_1_level = max(min(spoke_options.dim_1_level, 1), 0);
dim_2_level = max(min(spoke_options.dim_2_level, 1), 0);
dim_3_level = max(min(spoke_options.dim_3_level, 1), 0);
dim_4_level = max(min(spoke_options.dim_4_level, 1), 0);
axis_size = max(min(spoke_options.axis_size, 1), 0);
axis_color = spoke_options.axis_color;

dim_1_label = spoke_options.labels{1};
dim_2_label = spoke_options.labels{2};
dim_3_label = spoke_options.labels{3};
dim_4_label = spoke_options.labels{4};

title = spoke_options.title;

% The screen size
screen_rect = Screen('Rect', screen_ptr);
screen_cent = [round(screen_rect(3)/2) round(screen_rect(4)/2)];

% Spoke --------------------------

% Axes
tmp_axis_length = round(axis_size*min(screen_rect(3:4)));

tmp_vert_axis_x = screen_cent(1);
tmp_vert_axis_start_y = round((screen_rect(4) - tmp_axis_length)/2);
tmp_vert_axis_end_y = tmp_vert_axis_start_y + tmp_axis_length;
vert_axis = ...
    [tmp_vert_axis_x tmp_vert_axis_start_y ...
    tmp_vert_axis_x tmp_vert_axis_end_y];

Screen('DrawLine', screen_ptr, axis_color, ...
    vert_axis(1), vert_axis(2), vert_axis(3), vert_axis(4), pen_width);

tmp_horiz_axis_y = screen_cent(2);
tmp_horiz_axis_start_x = round((screen_rect(3) - tmp_axis_length)/2);
tmp_horiz_axis_end_x = tmp_horiz_axis_start_x + tmp_axis_length;
horiz_axis = ...
    [tmp_horiz_axis_start_x tmp_horiz_axis_y ...
    tmp_horiz_axis_end_x tmp_horiz_axis_y];

Screen('DrawLine', screen_ptr, axis_color, ...
    horiz_axis(1), horiz_axis(2), horiz_axis(3), horiz_axis(4), pen_width);

% Ticks
tick_length = default_tick_size_ratio*min(screen_rect(3:4));

for t = 0:default_num_ticks-1

    Screen('DrawLine', screen_ptr, axis_color, ...
        round(vert_axis(1) - tick_length/2), ...
        round(vert_axis(2) + ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))/2), ...
        round(vert_axis(1) + tick_length/2), ...
        round(vert_axis(2) + ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))/2), ...
        pen_width);

    Screen('DrawLine', screen_ptr, axis_color, ...
        round(vert_axis(1) - tick_length/2), ...
        round(vert_axis(4) - ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))/2), ...
        round(vert_axis(1) + tick_length/2), ...
        round(vert_axis(4) - ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))/2), ...
        pen_width);

    Screen('DrawLine', screen_ptr, axis_color, ...
        round(horiz_axis(1) + ...
        t*(1/default_num_ticks)*(horiz_axis(3) - horiz_axis(1))/2), ...
        round(horiz_axis(2) - tick_length/2), ...
        round(horiz_axis(1) + ...
        t*(1/default_num_ticks)*(horiz_axis(3) - horiz_axis(1))/2), ...
        round(horiz_axis(2) + tick_length/2), ...
        pen_width);

    Screen('DrawLine', screen_ptr, axis_color, ...
        round(horiz_axis(3) - ...
        t*(1/default_num_ticks)*(horiz_axis(3) - horiz_axis(1))/2), ...
        round(horiz_axis(2) - tick_length/2), ...
        round(horiz_axis(3) - ...
        t*(1/default_num_ticks)*(horiz_axis(3) - horiz_axis(1))/2), ...
        round(horiz_axis(2) + tick_length/2), ...
        pen_width);

end

% The lines

% Dim 1 to 2
dim_12_line = round([vert_axis(1), ...
    (1 - dim_1_level)*screen_cent(2) + dim_1_level*vert_axis(2), ...
    (1 - dim_2_level)*screen_cent(1) + dim_2_level*horiz_axis(3), ...
    horiz_axis(2)]);
Screen('DrawLine', screen_ptr, axis_color, ...
    dim_12_line(1), dim_12_line(2), dim_12_line(3), dim_12_line(4), ...
    pen_width);

% Dim 2 to 3
dim_23_line = round([(1 - dim_2_level)*screen_cent(1) + dim_2_level*horiz_axis(3), ...
    horiz_axis(2), ...
    vert_axis(1), ...
    (1 - dim_3_level)*screen_cent(2) + dim_3_level*vert_axis(4)]);
Screen('DrawLine', screen_ptr, axis_color, ...
    dim_23_line(1), dim_23_line(2), dim_23_line(3), dim_23_line(4), ...
    pen_width);

% Dim 3 to 4
dim_34_line = round([vert_axis(1), ...
    (1 - dim_3_level)*screen_cent(2) + dim_3_level*vert_axis(4), ...
    (1 - dim_4_level)*screen_cent(1) + dim_4_level*horiz_axis(1), ...
    horiz_axis(2)]);
Screen('DrawLine', screen_ptr, axis_color, ...
    dim_34_line(1), dim_34_line(2), dim_34_line(3), dim_34_line(4), ...
    pen_width);

% Dim 4 to 1
dim_41_line = round([(1 - dim_4_level)*screen_cent(1) + dim_4_level*horiz_axis(1), ...
    horiz_axis(2), ...
    vert_axis(1), ...
    (1 - dim_1_level)*screen_cent(2) + dim_1_level*vert_axis(2)]);
Screen('DrawLine', screen_ptr, axis_color, ...
    dim_41_line(1), dim_41_line(2), dim_41_line(3), dim_41_line(4), ...
    pen_width);
% Text ---------------------------

% The labels
[old_font_name old_font_number] = ...
    Screen('TextFont', screen_ptr, default_font);
old_font_size = ...
    Screen('TextSize', screen_ptr, default_font_size);

% Dim 1
[bounds tmp] = Screen('TextBounds', screen_ptr, dim_1_label, -100, -100);
Screen('DrawText', screen_ptr, dim_1_label, ...
    round(vert_axis(1) - bounds(3)/2), ...
    round(vert_axis(2) - bounds(4) - text_label_buffer*bounds(4)), ...
    axis_color);

% Dim 3
[bounds tmp] = Screen('TextBounds', screen_ptr, dim_3_label, -100, -100);
Screen('DrawText', screen_ptr, dim_3_label, ...
    round(vert_axis(1) - bounds(3)/2), ...
    round(vert_axis(4) + text_label_buffer*bounds(4)), ...
    axis_color);

% Dim 2
[bounds tmp] = Screen('TextBounds', screen_ptr, dim_2_label, -100, -100);
Screen('DrawText', screen_ptr, dim_2_label, ...
    round(horiz_axis(3) + text_label_buffer*bounds(4)), ...
    round(round(horiz_axis(2) - bounds(4)/2)), ...
    axis_color);

% Dim 4
[bounds tmp] = Screen('TextBounds', screen_ptr, dim_4_label, -100, -100);
Screen('DrawText', screen_ptr, dim_4_label, ...
    round(horiz_axis(1) - bounds(3) - text_label_buffer*bounds(4)), ...
    round(round(horiz_axis(2) - bounds(4)/2)), ...
    axis_color);

% Title
Screen('DrawText', screen_ptr, title, ...
    round(horiz_axis(1) - bounds(3) - text_label_buffer*bounds(4)), ...
    round(vert_axis(2) - bounds(4) - text_label_buffer*bounds(4)), ...
    axis_color);

Screen('TextFont', screen_ptr, old_font_name);
Screen('TextSize', screen_ptr, old_font_size);

% Finalization  -----------------------------------------------

spoke_locs.dim_12_line = dim_12_line;
spoke_locs.dim_23_line = dim_23_line;
spoke_locs.dim_34_line = dim_34_line;
spoke_locs.dim_41_line = dim_41_line;
spoke_locs.vert_axis = vert_axis;
spoke_locs.horiz_axis = horiz_axis;
