function car_texture = RDCL_DrawCar(angles, locs)

% function RDCL_DrawCar
%
% Draws a simple car with 2 dimensions, angle of back and front of cab.
%
% Arguments:
%    angles = 1x2 vector [front_angle back_angle]. In deg. 0 is straight
%      0-65 degrees.
%    up. Positive moves the angle toward the middle of the car.
%      front_angle: angle of front of cab
%      back_angle: angle of back of cab
%    locs = 1x2 vector [x_loc y_loc] (optional)
%      for upper left corner of car
%      if not given, will default to 0, 0
%    screen_ptr = screen pointer
%
% Returns:
%    car_texture = a texture of the car.
%
% Example:
%   car_texture = RDCL_DrawCar([45 45], [], screen_ptr);
%   Screen('DrawTexture', screen_ptr, car_texture)
%   Screen('Flip', screen_ptr);
%
% Andrew Cohen
% 12/10/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% 5/20/19 alc small changes. checked that it still works.
% 5/24/19 alc made screen_ptr a global

global screen_ptr

% Error checking
if angles(1) < 0 || angles(1) > 65 || angles(2) < 0 || angles(2) > 65
    error('RDCL_DrawCar: Angles must be between 0 and 65 deg');
end

% Constants
car_body_color = [58 137 216];
window_color = [132 132 132];
outline_color = [56 70 83];

% Arguments
front_angle = angles(1)*pi/180; % rads
back_angle = angles(2)*pi/180; % rads

if length(locs == 2)
    x_loc = locs(1);
    y_loc = locs(2);
else
    x_loc = 0;
    y_loc = 0;
end

% Read in the car base
car_base = imread('./RDCL_Functions/RDCL_Resources/car_base.png');

% Get size
image_width = size(car_base, 2);
image_height = size(car_base, 1);

% Make the texture an offscreen buffer so can draw other things into it
car_texture = Screen('OpenOffscreenWindow', screen_ptr);

% Put the car base into the texture
Screen('PutImage', car_texture, car_base, ...
    [x_loc, y_loc, x_loc + image_width, y_loc + image_height]);

% Locations of corners
back_roof_corner_x = x_loc + 48 + 37/tan(pi/2 - back_angle);
back_roof_corner_y = y_loc + 1;

front_roof_corner_x = x_loc + 270 - 49/tan(pi/2 - front_angle);
front_roof_corner_y = y_loc + 1;

back_window_corner_x = x_loc + 73 + 35/tan(pi/2 - back_angle);
back_window_corner_y = y_loc + 6;

front_window_corner_x = x_loc + 255 - 42/tan(pi/2 - front_angle);
front_window_corner_y = y_loc + 6;

% Draw the roof
roof_rect = [...
    x_loc + 48 y_loc + 38; ...
    back_roof_corner_x back_roof_corner_y; ...
    front_roof_corner_x front_roof_corner_y; ...
    x_loc + 270 y_loc + 50; ...
    x_loc + 255 y_loc + 48; ...
    front_window_corner_x front_window_corner_y; ...
    back_window_corner_x back_window_corner_y; ...
    x_loc + 73 y_loc + 44 ];

Screen('FillPoly', car_texture, car_body_color, roof_rect);

Screen('DrawLine', car_texture, outline_color, x_loc + 48, y_loc + 38, ...
    back_roof_corner_x, back_roof_corner_y);
Screen('DrawLine', car_texture, outline_color, back_roof_corner_x, ...
    back_roof_corner_y, front_roof_corner_x, front_roof_corner_y);
Screen('DrawLine', car_texture, outline_color, front_roof_corner_x, ...
    front_roof_corner_y, x_loc + 270, y_loc + 50);

% Draw the windows
back_window_rect = [...
    x_loc + 73 y_loc + 44; ...
    back_window_corner_x back_window_corner_y; ...
    x_loc + 155 back_window_corner_y; ...
    x_loc + 158 y_loc + 44; ...
    ];

Screen('FillPoly', car_texture, window_color, back_window_rect);

front_window_rect = [...
    x_loc + 166 front_window_corner_y; ...
    front_window_corner_x front_window_corner_y; ...
    x_loc + 255 y_loc + 48; ...
    x_loc + 171 y_loc + 44; ...
    ];

Screen('FillPoly', car_texture, window_color, front_window_rect);

Screen('DrawLine', car_texture, outline_color, x_loc + 73, y_loc + 44, ...
    back_window_corner_x, back_window_corner_y);
Screen('DrawLine', car_texture, outline_color, back_window_corner_x, ...
    back_window_corner_y, front_window_corner_x, front_window_corner_y);
Screen('DrawLine', car_texture, outline_color, front_window_corner_x, ...
    front_window_corner_y, x_loc + 255, y_loc + 49);
