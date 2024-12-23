
function RDCL_SaveFigure(file_name)

% function RDCL_SaveFigure
%
% Saves the current figure to a pdf file.
%
% Arguments:
%    file_name = the name of the saved file (including path). The name
%      should not include the file type extension.
%
% Returns:
%    none
%
% Example:
%   RDCL_SaveFigure('tmp');
%
% Andrew Cohen
% 4/6/09
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%   Can be called from the command line.
%   Does not allow differnet file types yet.
%
% Known Bugs:
%   none

% Pixels per inch
ppi = 72;

% Change the position mode to automatic, so things adjust correctly
set(gcf, 'PaperPositionMode', 'auto');

% Change the paper size to fit the figure
position = get(gcf, 'Position');
set(gcf, 'PaperSize', [position(3)/ppi position(4)/ppi]);

% Print it
print('-dpdf', '-r300', file_name);
% print('-depsc', '-tiff', '-r300', file_name);
