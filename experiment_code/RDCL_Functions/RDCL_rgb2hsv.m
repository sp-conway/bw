function hsv = RDCL_rgb2hsv(rgb)

%
% Function: RDCL_rgb2hsb
%
% Converts red/green/blue color to hue/saturation/brightness(value) color.
%
% Arguments:
%   rgb = 1x3 vector of [red green blue] (0-255)
%
% Return:
%   hsv = 1x3 vector or [hue (0-360) saturation (0-100) value (0-100)]
%
% Example:
%   hsv = RDCL_rgb2hsv([125 0 255]);
%
% Notes:
%   From http://delphi.about.com/od/adptips2006/qt/RgbToHsb.htm
%
% Andrew Cohen
% 6/25/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants
r = rgb(1);
g = rgb(2);
b = rgb(3);

% Convert
minRGB = min(min(r, g), b);
maxRGB = max(max(r, g), b);

delta = maxRGB - minRGB;

v = maxRGB;

if maxRGB ~= 0
    s = 255 * delta / maxRGB;
else
    s = 0;
end

if s ~= 0
    if r == maxRGB
        h = (g - b)/delta;
    elseif g == maxRGB
        h = 2 + (b - r)/delta;
    elseif b == maxRGB
        h = 4 + (r - g)/delta;
    end
else
    h = -1;
end

h = h * 60;

if h < 0
    h = h + 360;
end

s = s*100/255;
v = v*100/255;

hsv = [h s v];
