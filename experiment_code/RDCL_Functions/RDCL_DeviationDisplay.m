function bar_locs = RDCL_DeviationDisplay(deviation_options)

% Displays a deviation graph
%
%  screen_ptr = a pointer to the psychtoolbox screen
%   deviation_options = a struct with the arguments to make the bars
%     .dim_1_level: 0 to 1
%     .dim_2_level: 0 to 1
%     .dim_3_level: 0 to 1
%     .dim_4_level: 0 to 1
%     .axis_size: 0 = smallest, 1 = takes up whole screen
%         Uses the smaller of the x and y dimensions for max size
%     .axis_color: [r g b] color
%     .labels: cell of 4 labels (left to right)
%     .title: overall title
%
%  bar_rects = a struct with the locations of the bars
%    .dim_1_bar
%    .dim_2_bar
%    .dim_3_bar
%    .dim_4_bar
%    .vert_axis
%    .horiz_axis
%
% Example
%
%  deviation_options.dim_1_level = .75;
%  deviation_options.dim_2_level = .25;
%  deviation_options.dim_3_level = .15;
%  deviation_options.dim_4_level = .95;
%  deviation_options.axis_size = .7;
%  deviation_options.axis_color = [255 255 255];
%  deviation_options.labels = {'AGE', 'TOT. MILES', '10K', 'MOT.'};
%  deviation_options.title = 'Runner = 1';
%
%  RDCL_DeviationDisplay(deviation_options)
%
% Notes:
%   dim_1 is left, moving right
%
% Andrew Cohen
% 1/14/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change log:
%   5/24/19 alc added to RDCL functions.

% Constants -------------------------------------------

global screen_ptr;

% Defaults
default_font = 'Arial';
default_font_size = 26;

default_num_ticks = 10; % Not incl. origin
default_tick_size_ratio = .01;

default_axis_w_h_ratio = 1.5;

default_bar_width_ratio = .1; % Prop of horiz axis

pen_width = 1;

text_label_buffer = 1.5; % of vert extent

% Get options
dim_1_level = max(min(deviation_options.dim_1_level, 1), 0);
dim_2_level = max(min(deviation_options.dim_2_level, 1), 0);
dim_3_level = max(min(deviation_options.dim_3_level, 1), 0);
dim_4_level = max(min(deviation_options.dim_4_level, 1), 0);

dim_1_level = 2*dim_1_level - 1; % Change to -1 to 1
dim_2_level = 2*dim_2_level - 1;
dim_3_level = 2*dim_3_level - 1;
dim_4_level = 2*dim_4_level - 1;

axis_size = max(min(deviation_options.axis_size, 1), 0);
axis_color = deviation_options.axis_color;

dim_1_label = deviation_options.labels{1};
dim_2_label = deviation_options.labels{2};
dim_3_label = deviation_options.labels{3};
dim_4_label = deviation_options.labels{4};

title = deviation_options.title;

% The screen size
screen_rect = Screen('Rect', screen_ptr);
screen_cent = [round(screen_rect(3)/2) round(screen_rect(4)/2)];

% Bars --------------------------

% Axes
tmp_axis_h_length = round(axis_size*min(screen_rect(3:4)));
tmp_axis_v_length = round((1/default_axis_w_h_ratio)*tmp_axis_h_length);

tmp_vert_axis_x = round(screen_cent(1) - tmp_axis_h_length/2);
tmp_vert_axis_start_y = round((screen_rect(4) - tmp_axis_v_length)/2);
tmp_vert_axis_end_y = tmp_vert_axis_start_y + tmp_axis_v_length;
vert_axis = ...
    [tmp_vert_axis_x tmp_vert_axis_start_y ...
    tmp_vert_axis_x tmp_vert_axis_end_y];

Screen('DrawLine', screen_ptr, axis_color, ...
    vert_axis(1), vert_axis(2), vert_axis(3), vert_axis(4), pen_width);

tmp_horiz_axis_y = screen_cent(2);
tmp_horiz_axis_start_x = round(screen_cent(1) - tmp_axis_h_length/2);
tmp_horiz_axis_end_x = tmp_horiz_axis_start_x + tmp_axis_h_length;
horiz_axis = ...
    [tmp_horiz_axis_start_x tmp_horiz_axis_y ...
    tmp_horiz_axis_end_x tmp_horiz_axis_y];

Screen('DrawLine', screen_ptr, axis_color, ...
    horiz_axis(1), horiz_axis(2), horiz_axis(3), horiz_axis(4), pen_width);

% Ticks
tick_length = default_tick_size_ratio*min(screen_rect(3:4));

for t = 0:default_num_ticks

    Screen('DrawLine', screen_ptr, axis_color, ...
        round(vert_axis(1) - tick_length/2), ...
        round(vert_axis(2) + ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))), ...
        round(vert_axis(1) + tick_length/2), ...
        round(vert_axis(2) + ...
        t*(1/default_num_ticks)*(vert_axis(4) - vert_axis(2))), ...
        pen_width);

end

% The Bars

bar_width = round(default_bar_width_ratio*tmp_axis_h_length);

for d = 1:4

    switch(d)
        case 1, dim_level = dim_1_level;
        case 2, dim_level = dim_2_level;
        case 3, dim_level = dim_3_level;
        case 4, dim_level = dim_4_level;
    end

    bar_x_start = round(horiz_axis(1) + d*tmp_axis_h_length/5 - bar_width/2);
    bar_x_end = bar_x_start + bar_width;

    if dim_level >= 0

        bar_y_start = round(horiz_axis(2) - dim_level*tmp_axis_v_length/2);
        bar_y_end = horiz_axis(2);

    elseif dim_level < 0

        bar_y_start = horiz_axis(2);
        bar_y_end = round(horiz_axis(2) - dim_level*tmp_axis_v_length/2);

    end

    bar = [bar_x_start bar_y_start bar_x_end bar_y_end];

    Screen('FillRect', screen_ptr, axis_color, bar);

    switch(d)
        case 1, dim_1_bar = bar;
        case 2, dim_2_bar = bar;
        case 3, dim_3_bar = bar;
        case 4, dim_4_bar = bar;
    end

end

% Text ---------------------------

% The labels
[old_font_name old_font_number] = ...
    Screen('TextFont', screen_ptr, default_font);
old_font_size = ...
    Screen('TextSize', screen_ptr, default_font_size);

for d = 1:4

    switch(d)
        case 1
            curr_bar = dim_1_bar;
            curr_label = dim_1_label;
        case 2
            curr_bar = dim_2_bar;
            curr_label = dim_2_label;
        case 3
            curr_bar = dim_3_bar;
            curr_label = dim_3_label;
        case 4
            curr_bar = dim_4_bar;
            curr_label = dim_4_label;
    end

    [bounds tmp] = Screen('TextBounds', screen_ptr, curr_label, -100, -100);
    Screen('DrawText', screen_ptr, curr_label, ...
        round(horiz_axis(1) + d*tmp_axis_h_length/5 - bounds(3)/2), ...
        round(vert_axis(4) + text_label_buffer*bounds(4)), ...
        axis_color);

end

% Title
[bounds tmp] = Screen('TextBounds', screen_ptr, title, -100, -100);
Screen('DrawText', screen_ptr, title, ...
    round(screen_cent(1) - bounds(3)/2), ...
    round(vert_axis(4) + bounds(4) + 3*text_label_buffer*bounds(4)), ...
    axis_color);

Screen('TextFont', screen_ptr, old_font_name);
Screen('TextSize', screen_ptr, old_font_size);

% Finalization  -----------------------------------------------

bar_locs.dim_1_bar = dim_1_bar;
bar_locs.dim_2_bar = dim_2_bar;
bar_locs.dim_3_bar = dim_3_bar;
bar_locs.dim_4_bar = dim_4_bar;
bar_locs.vert_axis = vert_axis;
bar_locs.horiz_axis = horiz_axis;
