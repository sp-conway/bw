function perturbedDotPattern = RDCL_PerturbDotPattern(dotPattern, params, varargin)

%
% Function: RDCL_PerturbDotPattern
%
% Perturbs a dot pattern, i.e., generates a random dot pattern from a
%   prototype.
%
% Arguments:
%   dotPattern = the prototype to perturb.
%   params = a structure containing the type and amount of perturbation.
%     params{1} = 'Unif' or 'Norm'. Uniform or normal distribution.
%     params{2} = If 'Unif' = [min_horiz max_horiz min_vert max_vert] or
%                   [max_horiz max_vert]
%                   (in either direction, so basically doubles range)
%               = If 'Norm' = [mu_horiz sigma_horiz min_horiz max_horiz ...
%                     mu_vert sigma_vert min_vert max_vert] or
%                   [mu_horiz sigma_horiz max_horiz ...
%                     mu_vert sigma_vert max_vert] or
%                   [mu_horiz sigma_horiz mu_vert sigma_vert] or
%                   [sigma_horiz sigma_vert]
%                   (in pixels, covers both directions, centered at point).
%   varargin:
%     'Clip' = clip at orginal grid.
%
% Return:
%   perturbedDotPattern = a structure containing
%     perturbedDotPattern{1} = a matrix of the dot positions in pixels.
%       The lower-left corner of the grid is (0, 0).
%     perturbedDotPattern{2} = the gridSize.
%
% Example:
%   dotPattern = RDCL_GenerateDotPattern([100 100], 9);
%   perturbedDotPattern = ...
%     RDCL_PerturbDotPattern(dotPattern, {'Unif', [5 20 5 20]}, 'Clip');
%   perturbedDotPattern = ...
%     RDCL_PerturbDotPattern(dotPattern, {'Norm', [0 20 -Inf Inf 0 20 -Inf Inf]});
%
% Notes:
%   Assumes dimensions are independent.
%   mu in 'Norm' is mean from 0 in either direction.
%   gridSize is unchanged, even if points go outside grid.
%   'Norm', point can go on either side, min/max based on which side goes to.
%
% Andrew Cohen
% 6/24/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants
nDots = size(dotPattern{1}, 1);

% Defaults
clip = 0;

% Get user options
if nargin > 0
    for k = 1:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'Clip')
            clip = 1;
        else
            error('Incorrect RDCL_PerturbDotPattern option');
        end

    end
end

% Perturb the points
if strcmp(params{1}, 'Unif')

    if length(params{2}) == 4
        min_h = params{2}(1); max_h = params{2}(2);
        min_v = params{2}(3); max_v = params{2}(4);
    else
        min_h = 0; max_h = params{2}(1);
        min_v = 0; max_v = params{2}(2);
    end

    for p = 1:nDots

        for d = 1:2 % 2 dimensions

            % Horizontal
            while(1)

                % How much to perturb
                if d == 1
                    min_p = min_h; max_p = max_h;
                else
                    min_p = min_v; max_p = max_v;
                end

                r = round(min_p + (max_p - min_p)*rand);
                r = r * (2*(rand < .5) - 1); % sign

                % Do it
                perturbedDotPattern{1}(p, d) = dotPattern{1}(p, d) + r;

                % Clip
                if clip
                    if perturbedDotPattern{1}(p, d) >= 0 && perturbedDotPattern{1}(p, d) <= dotPattern{2}(d)
                        break;
                    end
                else
                    break;
                end

            end

        end

    end

elseif strcmp(params{1}, 'Norm')

    if length(params{2}) == 8
        mu_horiz = params{2}(1); sigma_horiz = params{2}(2);
        min_horiz = params{2}(3); max_horiz = params{2}(4);
        mu_vert = params{2}(5); sigma_vert = params{2}(6);
        min_vert = params{2}(7); max_vert = params{2}(8);
    elseif length(params{2}) == 6
        mu_horiz = params{2}(1); sigma_horiz = params{2}(2);
        min_horiz = -Inf; max_horiz = params{2}(3);
        mu_vert = params{2}(4); sigma_vert = params{2}(5);
        min_vert = -Inf; max_vert = params{2}(6);
    elseif length(params{2}) == 4
        mu_horiz = params{2}(1); sigma_horiz = params{2}(2);
        min_horiz = -Inf; max_horiz = Inf;
        mu_vert = params{2}(3); sigma_vert = params{2}(4);
        min_vert = -Inf; max_vert = Inf;
    elseif length(params{2}) == 2
        mu_horiz = 0; sigma_horiz = params{2}(1);
        min_horiz = -Inf; max_horiz = Inf;
        mu_vert = 0; sigma_vert = params{2}(2);
        min_vert = -Inf; max_vert = Inf;
    end

    for p = 1:nDots

        for d = 1:2 % 2 dimensions

            % Horizontal
            while(1)

                % How much to perturb
                if d == 1
                    mu = mu_horiz; sigma = sigma_horiz;
                    min_p = min_horiz; max_p = max_horiz;
                else
                    mu = mu_vert; sigma = sigma_vert;
                    min_p = min_vert; max_p = max_vert;
                end

                which_sign = (2*(rand < .5) - 1); % sign
                mu = mu * which_sign;

                while(1)
                    r = round(mu + sigma*randn);

                    if which_sign == 1
                        if r >= min_p && r <= max_p
                            break;
                        end
                    else
                        if r <= -min_p && r >= -max_p
                            break;
                        end
                    end
                end

                % Do it
                perturbedDotPattern{1}(p, d) = dotPattern{1}(p, d) + r;

                % Clip
                if clip
                    if perturbedDotPattern{1}(p, d) >= 0 && perturbedDotPattern{1}(p, d) <= dotPattern{2}(d)
                        break;
                    end
                else
                    break;
                end

            end

        end

    end

else

    error('Incorrect param type in RDCL_PerturbDotPattern');

end

% Keep grid size
perturbedDotPattern{2} = dotPattern{2};
