function RDCL_DrawArc(arc_rect, arc_type, color)

% Draws an arc
%
% arc_rect: the bounding rectangle
% arc_type: 0=full oval, 1=lower half, 2=upper half
% color = pen color
%
% Andrew Cohen
% 1/15/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change log:
%   5/24/19 alc makde screen_ptr a global.

global screen_ptr;

% Constants
delta_dep = 1;

pen_width = 1;

% Arc constants
a = round((arc_rect(3) - arc_rect(1))/2);
b = round((arc_rect(4) - arc_rect(2))/2);

x0 = round((arc_rect(3) + arc_rect(1))/2);
y0 = round((arc_rect(4) + arc_rect(2))/2);

% x as independent variable
indep = x0 - a:delta_dep:x0 + a;
cents = [x0 y0];
consts = [a b];

% Compute y
[dep indep] = indep_2_dep(indep, cents, consts, arc_type);

% Draw it
Screen('DrawDots', screen_ptr, round([indep; dep]), pen_width, color);

% y as independent variable
if arc_type == 0
    indep = y0 - b:delta_dep:y0 + b;
elseif arc_type == 1
    indep = y0:delta_dep:y0 + b;
else
    indep = y0 - b:delta_dep:y0;
end
cents = [y0 x0];
consts = [b a];

% Compute x
[dep indep] = indep_2_dep(indep, cents, consts, 0);

% Draw it
Screen('DrawDots', screen_ptr, round([dep; indep]), pen_width, color);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dep indep] = indep_2_dep(indep, cents, consts, arc_type)

dep = [];

if arc_type == 0 || arc_type == 1

    dep = ...
        [dep cents(2) + ...
        consts(2).*sqrt(1 - (indep - cents(1)).^2./consts(1)^2)];

end

if arc_type == 0 || arc_type == 2

    dep = ...
        [dep cents(2) - ...
        consts(2).*sqrt(1 - (indep - cents(1)).^2/consts(1).^2)];

end

if arc_type == 0

    indep = [indep indep];

end

end
