function rgb = RDCL_hsv2rgb(hsv)

%
% Function: RDCL_hsv2rgb
%
% Converts hue/saturation/brightness(value) to red/green/blue color color.
%
% Arguments:
%   hsv = 1x3 vector or [hue (0-360) saturation (0-100) value (0-100)]
%
% Return:
%   rgb = 1x3 vector of [red green blue] (0-255)
%
% Example:
%   hsv = RDCL_rgb2hsv([125 0 255]);
%
% Notes:
%   From Foley, et al, 1994
%
% Andrew Cohen
% 6/25/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants
h = hsv(1);
s = hsv(2)/100;
v = hsv(3)/100;

% Convert
if v == 0
    r = 0;
    g = 0;
    b = 0;
else
    if s == 0
        r = v;
        g = v;
        b = v;
    else
        if h == 360
            h = 0;
        end

        h = h / 60;
        i = floor(h);
        f = h - i;
        p = v*(1 - s);
        q = v*(1 - s*f);
        t = v*(1 - s*(1 - f));

        switch i
            case 0
                r = v; g = t; b = p;
            case 1
                r = q; g = v; b = p;
            case 2
                r = p; g = v; b = t;
            case 3
                r = p; g = q; b = v;
            case 4
                r = t; g = p; b = v;
            case 5
                r = v; g = p; b = q;
        end

    end
end

r = r*255;
g = g*255;
b = b*255;

rgb = [r g b];
