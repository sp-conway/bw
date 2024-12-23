function [ball_1_before_x ball_1_after_x ball_1_y ball_2_before_x ball_2_after_x ball_2_y mr ln_mr] = ...
    RDCL_TwoBallCollision(varargin)

%
% Function: RDCL_TwoBallCollision
%
% Shows a collision of 2 balls. Velocities, mass durations, sizes, and ball colors
%   can be vectors in which case there will be multiple collisions.
%
% Arguments:
%   varargin:
%     'Ball1BeforeX' = collison axis velocity Ball 1 before collision in
%       pixels/second.
%     'Ball1AfterX' = collison axis velocity Ball 1 after collision in
%       pixels/second.
%     'Ball2BeforeX' = collison axis velocity Ball 2 before collision in
%       pixels/second.
%     'Ball2AfterX' = collison axis velocity Ball 2 after collision in
%       pixels/second.
%     'Ball1Y' = sweep axis velocity Ball 1 before and after collision in
%       pixels/second.
%     'Ball2Y' = collison axis velocity Ball 2 before and after collision in
%       pixels/second.%
%     'MassRatio' = mass Ball 1 / mass Ball 2;
%     'LnMassRatio' = log(mass Ball 1 / mass Ball 2);
%     'DurationBefore' = duration of collision in seconds before collision;
%        Must be the same for all collisions.
%     'DurationAfter' = duration of collision in seconds after collision;
%        Must be the same for all collisions.
%     'Ball1Color' = Ball 1 color in RGB; Each row is a color.
%     'Ball2Color' = Ball 2 color in RGB; Each row is a color.
%     'bgColor' = background color in RGB;
%     'Ball1Size' = Ball 1 size in pixels;
%     'Ball2Size' = Ball 2 size in pixels;
%     'PauseAtStart' = How long to pause before the balls move.
%     'PauseAtEnd' = How long to pause when balls finish.
%     'TextDuring' = Text to display during collision at bottom of screen.
%     'TextAfter' = Text to display after collision  at bottom of screen..
%
% Return:
%   ball_1_before_x ball_1_after_x ball_1_y ball_2_before_x ball_2_after_x ball_2_y mr ln_mr =
%     the ball velocities and mr and log(mr).
%
% Example:
%   RDCL_TwoBallCollision('Ball1BeforeX', 100, 'Ball2BeforeX', -100,
%   'MassRatio', 2, 'Ball1Size', 100);
%
% Notes:
%   Balls 1 & 2 always collide in middle of the screen with Ball 1 on the
%     left and Ball 2 on the right.
%
% Andrew Cohen
% 5/11/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%
% 5/16/06: Now handles multiple collisions.

% Constants
global screen_ptr

[screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle();
[screen_width screen_height] = Screen('WindowSize', screen_ptr);
ifi = RDCL_GetIFI('debug', 'T');

gray = RDCL_GetGray();

% Default parameters
ball_1_before_x = NaN;
ball_1_after_x = NaN;
ball_1_y = 0;
ball_2_before_x = NaN;
ball_2_after_x = NaN;
ball_2_y = 0;
mr = NaN;
num_seconds_before = 1;
num_seconds_after = 1;
bg_color = [0 0 0];
ball_1_color = [120 0 0];
ball_2_color = [120 0 0];
ball_1_size = 75;
ball_2_size = 75;
pause_start_time = 0;
pause_end_time = 0;
text_to_show_during = '';

% Get user options
if nargin > 0
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'Ball1BeforeX')
            ball_1_before_x = varargin{k+1};
        elseif strcmp(option_text,'Ball1AfterX')
            ball_1_after_x = varargin{k+1};
        elseif strcmp(option_text, 'Ball2BeforeX')
            ball_2_before_x = varargin{k+1};
        elseif strcmp(option_text, 'Ball2AfterX')
            ball_2_after_x = varargin{k+1};
        elseif strcmp(option_text, 'Ball1Y')
            ball_1_y = varargin{k+1};
        elseif strcmp(option_text, 'Ball2Y')
            ball_2_y = varargin{k+1};
        elseif strcmp(option_text, 'MassRatio')
            mr = varargin{k+1};
        elseif strcmp(option_text, 'LnMassRatio')
            mr = exp(varargin{k+1});
        elseif strcmp(option_text, 'DurationBefore')
            num_seconds_before = varargin{k+1};
        elseif strcmp(option_text, 'DurationAfter')
            num_seconds_after = varargin{k+1};
        elseif strcmp(option_text, 'Ball1Color')
            ball_1_color = varargin{k+1};
        elseif strcmp(option_text, 'Ball2Color')
            ball_2_color = varargin{k+1};
        elseif strcmp(option_text, 'bgColor')
            bg_color = varargin{k+1};
        elseif strcmp(option_text, 'Ball1Size')
            ball_1_size = round(varargin{k+1});
        elseif strcmp(option_text, 'Ball2Size')
            ball_2_size = round(varargin{k+1});
        elseif strcmp(option_text, 'PauseAtStart')
            pause_start_time = varargin{k+1};
        elseif strcmp(option_text, 'PauseAtEnd')
            pause_end_time = varargin{k+1};
        elseif strcmp(option_text, 'TextDuring')
            text_to_show_during = varargin{k+1};
        else
            error('Incorrect display_collision option');
        end

    end
end

% Number of collisions
num_collisions = max([length(ball_1_before_x) length(ball_1_after_x) ...
    length(ball_2_before_x) length(ball_2_after_x)]);

% Error Checking
if ~isnan(ball_1_before_x) && ~isnan(ball_1_after_x) && ...
        ~isnan(ball_2_before_x) && ~isnan(ball_2_after_x) && ~isnan(mr)
    error('Cannot specify all 4 collision axis velocities and mass ratio');
elseif ~((~isnan(ball_1_before_x) && ~isnan(ball_2_before_x) && isnan(ball_1_after_x) && isnan(ball_2_after_x)) | ...
        (isnan(ball_1_before_x) && isnan(ball_2_before_x) && ~isnan(ball_1_after_x) && ~isnan(ball_2_after_x))) ...
        && ~isnan(mr)
    error('Need to specify exactly 2 before or after velocities when define mass ratio');
elseif mr <= 0
    error('Mass ratio must be positive');
elseif ball_1_size < 1 | ball_2_size < 1
    error('Ball size must be >= 1');
elseif num_seconds_before <= 0 | num_seconds_after <= 0
    error('Duration must be positive');
elseif length(num_seconds_before) > 1 | length(num_seconds_after) > 1
    error('Duration must be a scalar');
elseif size(bg_color, 1) > 1
    error('BG color must be a single row');
elseif size(ball_1_color, 1) ~= num_collisions && size(ball_1_color, 1) ~= 1
    error('Must specify 1 or the number of collisions for ball color');
elseif size(ball_2_color, 1) ~= num_collisions && size(ball_2_color, 1) ~= 1
    error('Must specify 1 or the number of collisions for ball color');
elseif size(ball_1_size, 1) ~= num_collisions && size(ball_1_size, 1) ~= 1
    error('Must specify 1 or the number of collisions for ball size');
elseif size(ball_2_size, 1) ~= num_collisions && size(ball_2_size, 1) ~= 1
    error('Must specify 1 or the number of collisions for ball size');
elseif ~isnan(ball_1_before_x) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(ball_2_before_x) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(ball_1_after_x) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(ball_2_after_x) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(ball_1_y) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(ball_2_y) && length(ball_1_before_x) ~= num_collisions
    error('All velocity vectors must have the same length');
elseif ~isnan(mr) && length(mr) ~= num_collisions
    error('Mass ratio must have the same length as velocities');
end

% Set 2 velocities given 2 and mass ratio
if ~isnan(mr) && ~isnan(ball_1_before_x) && ~isnan(ball_2_before_x)
    % If have entry velocities and mr

    ball_1_mass = ones(1, num_collisions);
    ball_2_mass = 1./mr;

	ball_1_after_x = (ball_1_mass - ball_2_mass)./(ball_1_mass + ball_2_mass).*...
		ball_1_before_x + 2*ball_2_mass./(ball_1_mass + ball_2_mass).*...
		ball_2_before_x;

	ball_2_after_x = 2*ball_1_mass./(ball_1_mass + ball_2_mass).*ball_1_before_x + ...
		(ball_2_mass - ball_1_mass)./(ball_1_mass + ball_2_mass).*ball_2_before_x;

elseif ~isnan(mr) && ~isnan(ball_1_after_x) & ~isnan(ball_2_after_x)
    % If have exit velocities and mr

    ball_1_mass = ones(1, num_collisions);
    ball_2_mass = 1./mr;

	ball_1_before_x = (ball_1_mass - ball_2_mass)./(ball_1_mass + ball_2_mass).*...
		ball_1_after_x + 2*ball_2_mass./(ball_1_mass + ball_2_mass).*...
		ball_2_after_x;

	ball_2_before_x = 2*ball_1_mass./(ball_1_mass + ball_2_mass).*ball_1_after_x + ...
		(ball_2_mass - ball_1_mass)./(ball_1_mass + ball_2_mass).*ball_2_after_x;

end

% Constants
frames_per_half_display_before = RDCL_Sec2Frames(num_seconds_before, ifi.duration);
frames_per_half_display_after = RDCL_Sec2Frames(num_seconds_after, ifi.duration);

% Reshape colors and sizes
num_to_repeat = num_collisions - size(ball_1_color, 1) + 1;
ball_1_color = repmat(ball_1_color, num_to_repeat, 1);
num_to_repeat = num_collisions - size(ball_2_color, 1) + 1;
ball_2_color = repmat(ball_2_color, num_to_repeat, 1);

num_to_repeat = num_collisions - length(ball_1_size) + 1;
ball_1_size = repmat(ball_1_size, 1, num_to_repeat);
num_to_repeat = num_collisions - length(ball_2_size) + 1;
ball_2_size = repmat(ball_2_size, 1, num_to_repeat);

% Figure out where balls collide (upper left corner)
collision_1_x = repmat(screen_middle_x, 1, num_collisions) - ball_1_size;
collision_1_y = [1:num_collisions].*repmat(screen_height/(num_collisions+1), 1, num_collisions) - round(ball_1_size./2);
collision_2_x = repmat(screen_middle_x, 1, num_collisions);
collision_2_y = [1:num_collisions].*repmat(screen_height/(num_collisions+1), 1, num_collisions) - round(ball_2_size./2);

% Figure out start and ending locations (upper left corner)
ball_1_start_x = collision_1_x - ball_1_before_x.*num_seconds_before;
ball_1_end_x = collision_1_x + ball_1_after_x.*num_seconds_after;
ball_1_start_y = collision_1_y - ball_1_y.*num_seconds_before;
ball_1_end_y = collision_1_y + ball_1_y.*num_seconds_after;

ball_2_start_x = collision_2_x - ball_2_before_x.*num_seconds_before;
ball_2_end_x = collision_2_x + ball_2_after_x.*num_seconds_after;
ball_2_start_y = collision_2_y - ball_2_y.*num_seconds_before;
ball_2_end_y = collision_2_y + ball_2_y.*num_seconds_after;

% Before collision
for t = 1:frames_per_half_display_before
    Screen('FillRect', screen_ptr, gray);

    for i = 1:num_collisions

        % Location
        ball_1_x_loc = ball_1_start_x(i) + t*(collision_1_x(i) - ball_1_start_x(i))*1/frames_per_half_display_before;
        ball_2_x_loc = ball_2_start_x(i) + t*(collision_2_x(i) - ball_2_start_x(i))*1/frames_per_half_display_before;
        ball_1_y_loc = ball_1_start_y(i) + t*(collision_1_y(i) - ball_1_start_y(i))*1/frames_per_half_display_before;
        ball_2_y_loc = ball_2_start_y(i) + t*(collision_2_y(i) - ball_2_start_y(i))*1/frames_per_half_display_before;

        Screen('FillOval', screen_ptr, ball_1_color(i, :), [ball_1_x_loc ball_1_y_loc ball_1_x_loc + ball_1_size(i) ball_1_y_loc + ball_1_size(i)]); % Ball 1
        Screen('FillOval', screen_ptr, ball_2_color(i, :), [ball_2_x_loc ball_2_y_loc ball_2_x_loc + ball_2_size(i) ball_2_y_loc + ball_2_size(i)]); % Ball 2

    end

    if ~strcmp('', text_to_show_during)
        [norm_bounds_rect offset_bounds_rect] = Screen('TextBounds', screen_ptr, text_to_show_during);
        Screen('DrawText', screen_ptr, text_to_show_during, screen_middle_x - round(norm_bounds_rect(3)/2), screen_height - norm_bounds_rect(4)*2);
    end

    Screen('Flip', screen_ptr);

    if t == 1
        WaitSecs(pause_start_time);
    end
end

% After collision
for t = 1:frames_per_half_display_after
    Screen('FillRect', screen_ptr, gray);

    for i = 1:num_collisions

        % Location
        ball_1_x_loc = collision_1_x(i) - t*(collision_1_x(i) - ball_1_end_x(i))*1/frames_per_half_display_after;
        ball_2_x_loc = collision_2_x(i) - t*(collision_2_x(i) - ball_2_end_x(i))*1/frames_per_half_display_after;
        ball_1_y_loc = collision_1_y(i) - t*(collision_1_y(i) - ball_1_end_y(i))*1/frames_per_half_display_after;
        ball_2_y_loc = collision_2_y(i) - t*(collision_2_y(i) - ball_2_end_y(i))*1/frames_per_half_display_after;

        Screen('FillOval', screen_ptr, ball_1_color(i, :), [ball_1_x_loc ball_1_y_loc ball_1_x_loc + ball_1_size(i) ball_1_y_loc + ball_1_size(i)]); % Ball 1
        Screen('FillOval', screen_ptr, ball_2_color(i, :), [ball_2_x_loc ball_2_y_loc ball_2_x_loc + ball_2_size(i) ball_2_y_loc + ball_2_size(i)]); % Ball 2

    end

    if ~strcmp('', text_to_show_during)
        [norm_bounds_rect offset_bounds_rect] = Screen('TextBounds', screen_ptr, text_to_show_during);
        Screen('DrawText', screen_ptr, text_to_show_during, screen_middle_x - round(norm_bounds_rect(3)/2), screen_height - norm_bounds_rect(4)*2);
    end

    Screen('Flip', screen_ptr);
end

WaitSecs(pause_end_time);

Screen('FillRect', screen_ptr, gray);
Screen('Flip', screen_ptr);

ln_mr = log(mr);
