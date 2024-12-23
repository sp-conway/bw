function RDCL_DrawProtractor(radius, angle, arc_color, line_color, fill, line_length, locs)

% function RDCL_DrawProtractor
%
% Draws a protractor-like stimulus.
%
%
% Arguments:
%    radius = scalar in pixels. Radius of the arc.
%    angle = scalar in degrees. Angle of the line. 0 = straight up.
%    arc_color = [rgb] (optional). Color of the arc lines. Default = white.
%    line_color = [rgb] (optional). Color of the lines or fill arc. Default = white.
%    fill = boolean (optional). A line or a filled arc. Default = false.
%    line_length = scalar in pixels (optional). The line/fill can have its own size.
%      Default = radius.
%    locs = 1x2 vector [x_loc y_loc] (optional)
%      for center of the protractor.
%      If not given, will default to 0, 0.
%      Assumes protractor is top half of circle.
%    screen_ptr = screen pointer
%
% Returns:
%    None.
%
% Andrew Cohen
% 5/4/18
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change log:
%   5/24/19 alc made screen_ptr a global

global screen_ptr

% Defaults
if (length(arc_color)==0) arc_color = [255 255 255]; end
if (length(line_color)==0) line_color = [255 255 255]; end
if (length(line_length)==0) line_length = radius; end
if (length(fill)==0) fill = false; end
if (length(locs)==0) locs = [0 0]; end

radius = round(radius);
line_length = round(line_length);

% The arc
startAngle = -90;
arcAngle = 180;
arcRect = [locs(1)-radius locs(2)-.5*radius locs(1)+radius locs(2)+1.5*radius];
Screen('DrawArc', screen_ptr, arc_color, arcRect, startAngle, arcAngle);
Screen('DrawLine', screen_ptr, arc_color, locs(1)-radius, locs(2)+.5*radius, locs(1)+radius, locs(2)+.5*radius, 1);

% The line
if (fill == true)
  diff = (line_length - radius)/2;
  lineRect = [locs(1)-line_length locs(2)-.5*line_length-diff locs(1)+line_length locs(2)+1.5*line_length-diff];
  Screen('FillArc', screen_ptr, line_color, lineRect, angle, 90-angle)
  Screen('DrawLine', screen_ptr, arc_color, locs(1), locs(2)+.5*radius, locs(1)+line_length, locs(2)+.5*radius, 1);
else
  diff = (line_length - radius)/2;
  startX = locs(1);
  startY = locs(2)+.5*line_length-diff;
  endX   = startX + line_length * sin(angle*pi/180);
  endY   = startY - line_length * cos(angle*pi/180) - diff;
  Screen('DrawLine', screen_ptr, line_color, startX, startY, endX, endY , 1);
end
