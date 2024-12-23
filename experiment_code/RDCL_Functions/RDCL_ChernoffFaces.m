function face_rects = RDCL_ChernoffFaces(loc, face_options)

% Displays a Chernoff face
%
%   loc = [x y] position of middle of face
%   face_options = a struct with the arguments to make the face
%     .smile_level: 0 = most frown, 1 = most smile
%     .mouth_height: 0 = most down, 1 = most up
%     .nose_size: 0 = smallest, 1 = largest
%     .eye_height: 0 = lowest, 1 = highest
%     .eye_width: 0 = closest, 1 = farthest
%     .pupil_pos: 0 = straight on, 1 = looking to right
%     .eyebrow_height: 0 = lowest, 1 = highest
%     .face_size: [x y] size of face rectangle
%     .face_color: [r g b] color of the face
%     .labels: cell of 4 labels (eye brows, eyes, nose, mouth)
%     .title: graph title
%     .show_text: 1 = yes, 0 = no
%  screen_ptr = a pointer to the psychtoolbox screen
%
%  face_rects = a struct with the rects for the face components
%    .mouth_rect
%    .nose_rect
%    .eye_X_rect (X is 1 or 2)
%    .eyebrow_y
%    .eyebrow_X_x (X is 1 or 2)
%    .face_rect
%
% Notes:
%
% Andrew Cohen
% 1/11/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change log:
%   5/24/19 alc made screen_ptr a global

% Constants -------------------------------------------

global screen_ptr

% Defaults
default_font = 'Arial';
default_font_size = 26;

text_x_spacing = 10;
text_line_spacing = 1.7;

pen_width = 1;

% Get options
smile_level = max(min(face_options.smile_level, 1), 0);
mouth_height = max(min(face_options.mouth_height, 1), 0);
nose_size = max(min(face_options.nose_size, 1), 0);
eye_height = max(min(face_options.eye_height, 1), 0);
eye_width = max(min(face_options.eye_width, 1), 0);
pupil_pos = max(min(face_options.pupil_pos, 1), 0);
eyebrow_height = max(min(face_options.eyebrow_height, 1), 0);
face_size = face_options.face_size;
face_color = face_options.face_color;

eyebrow_label = face_options.labels{1};
eye_label = face_options.labels{2};
nose_label = face_options.labels{3};
mouth_label = face_options.labels{4};

title = face_options.title;

show_text = face_options.show_text;

% Sizes

% Mouth
min_mouth_cent_y = .15;
max_mouth_cent_y = .3;
min_mouth_width_ratio = .4;
max_mouth_width_ratio = .5;

% Nose
min_nose_ratio = .15;
max_nose_ratio = .35;
nose_start_y_ratio = .33;
nose_width_2_length_ratio = .8;

% Eyes
eye_height_ratio = .15;
eye_width_ratio = .25;
min_eye_cent_y_ratio = .23;
max_eye_cent_y_ratio = nose_start_y_ratio;
min_eye_cent_x_ratio = .20;
max_eye_cent_x_ratio = .30;
eye_pupil_ratio = .65; % relative to eye height
eye_pupil_min_right = .25;
eye_pupil_max_right = .75;

% Eyebrows
eyebrow_width_ratio = eye_width_ratio;
eyebrow_min_y_ratio = .03; % proportion of face, but relative to eye position
eyebrow_max_y_ratio = .13;

% Figure out sizes ------------------------------------

% The screen size
pos_cent = loc;

% The face rect
tmp_face_length = face_size;
tmp_start_x = loc(1) - round(tmp_face_length)/2;
tmp_start_y = loc(2) - round(tmp_face_length)/2;

face_rect = [tmp_start_x tmp_start_y ...
    tmp_start_x + tmp_face_length tmp_start_y + tmp_face_length];

% The mouth rect
tmp_mouth_cent_x = pos_cent(1);
tmp_mouth_height =  min_mouth_cent_y + ...
    (max_mouth_cent_y - min_mouth_cent_y)*mouth_height;
tmp_mouth_cent_y = round(face_rect(4) - tmp_mouth_height*tmp_face_length);

tmp_mouth_abs_level = 2*abs(smile_level - .5);

tmp_mouth_width_ratio = min_mouth_width_ratio + ...
    (max_mouth_width_ratio - min_mouth_width_ratio)*...
    (1 - tmp_mouth_abs_level);
tmp_mouth_width = round(tmp_face_length*tmp_mouth_width_ratio);

tmp_max_mouth_width = round(max_mouth_width_ratio*tmp_face_length);

tmp_mouth_height = 2*(tmp_max_mouth_width - tmp_mouth_width);

tmp_mouth_start_x = round(tmp_mouth_cent_x - tmp_mouth_width/2);
tmp_mouth_end_x = round(tmp_mouth_cent_x + tmp_mouth_width/2);

tmp_mouth_start_y = round(tmp_mouth_cent_y - tmp_mouth_height/2);
tmp_mouth_end_y = round(tmp_mouth_cent_y + tmp_mouth_height/2);

mouth_rect = ...
    [tmp_mouth_start_x tmp_mouth_start_y tmp_mouth_end_x tmp_mouth_end_y];

if smile_level > .5
    mouth_rect(2) = round(mouth_rect(2) - tmp_mouth_height/4);
    mouth_rect(4) = round(mouth_rect(4) - tmp_mouth_height/4);
elseif smile_level < .5
    mouth_rect(2) = round(mouth_rect(2) + tmp_mouth_height/4);
    mouth_rect(4) = round(mouth_rect(4) + tmp_mouth_height/4);
end

% The nose rect
min_nose_length = round(tmp_face_length*min_nose_ratio);
max_nose_length = round(tmp_face_length*max_nose_ratio);

tmp_nose_length =  round(min_nose_length + ...
    (max_nose_length - min_nose_length)*nose_size);

tmp_nose_start_y = round(face_rect(2) + nose_start_y_ratio*tmp_face_length);
tmp_nose_end_y = round(tmp_nose_start_y + tmp_nose_length);

tmp_nose_width = round(tmp_nose_length*nose_width_2_length_ratio);
tmp_nose_start_x = round(pos_cent(1) - tmp_nose_width/2);
tmp_nose_end_x = round(pos_cent(1) + tmp_nose_width/2);

nose_rect = ...
    [tmp_nose_start_x tmp_nose_start_y tmp_nose_end_x tmp_nose_end_y];

% The eyes rect
tmp_eye_height = round(tmp_face_length*eye_height_ratio);
tmp_eye_width = round(tmp_face_length*eye_width_ratio);

tmp_eye_cent_y_ratio = ...
  min_eye_cent_y_ratio + (max_eye_cent_y_ratio - min_eye_cent_y_ratio)*(1-eye_height);
tmp_eye_start_y = round(face_rect(2) + ...
    tmp_face_length*tmp_eye_cent_y_ratio - tmp_eye_height/2);
tmp_eye_end_y = tmp_eye_start_y + tmp_eye_height;

eye_cent_x_ratio = min_eye_cent_x_ratio + ...
  (max_eye_cent_x_ratio - min_eye_cent_x_ratio)*(1-eye_width);

tmp_eye_1_start_x = round(face_rect(1) + ...
    tmp_face_length*eye_cent_x_ratio - tmp_eye_width/2);
tmp_eye_1_end_x = tmp_eye_1_start_x + tmp_eye_width;

eye_1_rect = [tmp_eye_1_start_x tmp_eye_start_y ...
    tmp_eye_1_end_x tmp_eye_end_y];

tmp_eye_2_start_x = round(face_rect(3) - ...
    tmp_face_length*eye_cent_x_ratio - tmp_eye_width/2);
tmp_eye_2_end_x = tmp_eye_2_start_x + tmp_eye_width;

eye_2_rect = [tmp_eye_2_start_x tmp_eye_start_y ...
    tmp_eye_2_end_x tmp_eye_end_y];

tmp_pupil_height = round(tmp_eye_height*eye_pupil_ratio);

tmp_pupil_start_y = round((tmp_eye_start_y + tmp_eye_end_y)/2 - ...
    tmp_pupil_height/2);
tmp_pupil_end_y = tmp_pupil_start_y + tmp_pupil_height;

tmp_pupil_x_ratio = (1 - pupil_pos)* eye_pupil_min_right + ...
    pupil_pos * eye_pupil_max_right;

tmp_pupil_1_start_x = round(tmp_pupil_x_ratio * ...
    (tmp_eye_1_end_x - tmp_eye_1_start_x) + tmp_eye_1_start_x - ...
    tmp_pupil_height/2);
tmp_pupil_1_end_x = tmp_pupil_1_start_x + tmp_pupil_height;

pupil_1_rect = [tmp_pupil_1_start_x tmp_pupil_start_y ...
    tmp_pupil_1_end_x tmp_pupil_end_y];

tmp_pupil_2_start_x = round(tmp_pupil_x_ratio * ...
    (tmp_eye_2_end_x - tmp_eye_2_start_x) + tmp_eye_2_start_x - ...
    tmp_pupil_height/2);
tmp_pupil_2_end_x = tmp_pupil_2_start_x + tmp_pupil_height;

pupil_2_rect = [tmp_pupil_2_start_x tmp_pupil_start_y ...
    tmp_pupil_2_end_x tmp_pupil_end_y];

% The eyebrows rect
tmp_eyebrow_width = round(eyebrow_width_ratio*tmp_face_length);

tmp_eyebrow_y_ratio = eyebrow_height*eyebrow_min_y_ratio + ...
    (1 - eyebrow_height)*eyebrow_max_y_ratio;
eyebrow_y = round(tmp_eyebrow_y_ratio*tmp_face_length + face_rect(2));

tmp_eyebrow_y_ratio = eyebrow_min_y_ratio + ...
  (eyebrow_max_y_ratio - eyebrow_min_y_ratio)*(1-eyebrow_height);
tmp_eyebrow_height = tmp_eyebrow_y_ratio*tmp_face_length;
eyebrow_y = eye_1_rect(2) - tmp_eyebrow_height;

tmp_eyebrow_1_start_x = round(tmp_eye_1_start_x + ...
    (tmp_eye_width - tmp_eyebrow_width)/2);
tmp_eyebrow_1_end_x = tmp_eyebrow_1_start_x + tmp_eyebrow_width;

eyebrow_1_x = [tmp_eyebrow_1_start_x tmp_eyebrow_1_end_x];

tmp_eyebrow_2_start_x = round(tmp_eye_2_start_x + ...
    (tmp_eye_width - tmp_eyebrow_width)/2);
tmp_eyebrow_2_end_x = tmp_eyebrow_2_start_x + tmp_eyebrow_width;

eyebrow_2_x = [tmp_eyebrow_2_start_x tmp_eyebrow_2_end_x];

% Draw the face ---------------------------------------

% The circle
%Screen('FrameArc', screen_ptr, face_color, face_rect, 0, 360, 1);
RDCL_DrawArc(face_rect, 0, face_color)

% The mouth
if smile_level > .5 && mouth_rect(2) ~= mouth_rect(4)
    %Screen('FrameArc', screen_ptr, face_color, mouth_rect, 90, 180, ...
    %    pen_width, pen_width);
    RDCL_DrawArc(mouth_rect, 1, face_color)
elseif smile_level < .5 && mouth_rect(2) ~= mouth_rect(4)
    %Screen('FrameArc', screen_ptr, face_color, mouth_rect, 270, 180, ...
    %    pen_width, pen_width);
    RDCL_DrawArc(mouth_rect, 2, face_color)
else
    Screen('DrawLine', screen_ptr, face_color, ...
        mouth_rect(1), mouth_rect(2), mouth_rect(3), mouth_rect(4), ...
        pen_width);
end

% The nose
Screen('DrawLine', screen_ptr, face_color, ...
    pos_cent(1), nose_rect(2), nose_rect(3), nose_rect(4), pen_width);
Screen('DrawLine', screen_ptr, face_color, ...
    nose_rect(3), nose_rect(4), nose_rect(1), nose_rect(4), pen_width);

% The eyes
% Screen('DrawArc', screen_ptr, face_color, eye_1_rect, 0, 360);
% Screen('DrawArc', screen_ptr, face_color, eye_2_rect, 0, 360);
%
% Screen('DrawArc', screen_ptr, face_color, pupil_1_rect, 0, 360);
% Screen('DrawArc', screen_ptr, face_color, pupil_2_rect, 0, 360);

RDCL_DrawArc(eye_1_rect, 0, face_color)
RDCL_DrawArc(eye_2_rect, 0, face_color)

RDCL_DrawArc(pupil_1_rect, 0, face_color)
RDCL_DrawArc(pupil_2_rect, 0, face_color)

% The eye brows
Screen('DrawLine', screen_ptr, face_color, ...
    eyebrow_1_x(1), eyebrow_y, eyebrow_1_x(2), eyebrow_y);
Screen('DrawLine', screen_ptr, face_color, ...
    eyebrow_2_x(1), eyebrow_y, eyebrow_2_x(2), eyebrow_y);

% Text --------------------------------------

if show_text == 1

    % Set fonts
    [old_font_name old_font_number] = ...
        Screen('TextFont', screen_ptr, default_font);
    old_font_size = ...
        Screen('TextSize', screen_ptr, default_font_size);

    % The labels
    [bounds tmp] = Screen('TextBounds', screen_ptr, 'Test', -100, -100);

    tmp_text_x = face_rect(3) + text_x_spacing;
    tmp_text_y = face_rect(2);

    tmp_label = ['Eyebrows = ' eyebrow_label];
    Screen('DrawText', screen_ptr, tmp_label, tmp_text_x, tmp_text_y, ...
        face_color);

    tmp_label = ['Eyes = ' eye_label];
    tmp_text_y = tmp_text_y + text_line_spacing*bounds(4);
    Screen('DrawText', screen_ptr, tmp_label, tmp_text_x, tmp_text_y, ...
        face_color);

    tmp_label = ['Nose = ' nose_label];
    tmp_text_y = tmp_text_y + text_line_spacing*bounds(4);
    Screen('DrawText', screen_ptr, tmp_label, tmp_text_x, tmp_text_y, ...
        face_color);

    tmp_label = ['Mouth = ' mouth_label];
    tmp_text_y = tmp_text_y + text_line_spacing*bounds(4);
    Screen('DrawText', screen_ptr, tmp_label, tmp_text_x, tmp_text_y, ...
        face_color);

    % Title
    Screen('DrawText', screen_ptr, title, tmp_text_x, ...
        face_rect(4) - bounds(4), face_color);

    % Reset fonts
    Screen('TextFont', screen_ptr, old_font_name);
    Screen('TextSize', screen_ptr, old_font_size);

end

% Finalization  -----------------------------------------------

face_rects.face_rect = face_rect;
face_rects.mouth_rect = mouth_rect;
face_rects.eye_1_rect = eye_1_rect;
face_rects.eye_2_rect = eye_2_rect;
face_rects.nose_rect = nose_rect;
face_rects.eyebrow_y = eyebrow_y;
face_rects.eyebrow_1_x = eyebrow_1_x;
face_rects.eyebrow_2_x = eyebrow_2_x;
