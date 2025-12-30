function color = branchColorVariant(baseColor, branchIdx)
%BRANCHCOLORVARIANT Lighten/darken the base color to distinguish branches.
%   COLOR = BRANCHCOLORVARIANT(BASECOLOR, BRANCHIDX) nudges BASECOLOR toward
%   white or black (depending on the branch index) so that multiple branches
%   that share the same sweep color remain visually distinct even when their
%   profiles overlap.

if nargin < 1 || isempty(baseColor)
    color = [0, 0, 0];
    return
end

color = double(baseColor(:).');
if numel(color) < 3
    color = repmat(color(1), 1, 3);
elseif numel(color) > 3
    color = color(1:3);
end
color = clampColor(color);

if nargin < 2 || ~isscalar(branchIdx) || ~isfinite(branchIdx) || branchIdx <= 1
    return
end

branchTier = max(0, branchIdx - 2);
strength = 0.25 + 0.15 * floor(branchTier / 2);
strength = min(max(strength, 0), 0.85);

if mod(branchIdx, 2) == 0
    % Even branches: lighten the color toward white.
    color = color + (1 - color) * strength;
else
    % Odd branches beyond the first: darken toward black.
    color = color * (1 - strength);
end
color = clampColor(color);
end

function value = clampColor(value)
value = max(0, min(1, value));
end
