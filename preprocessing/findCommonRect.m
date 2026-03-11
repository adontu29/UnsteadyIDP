function rect = findCommonRect(dataFolder, nErode)
% findCommonRect  Find an upstream rectangle that lies inside the valid PIV
% domain for ALL datasets in a folder.
%
% Rule (as requested):
%   1) Right edge at xR = xLE - 0.5*c  (conservative over all datasets)
%   2) Vertical bounds at zLEref +/- 0.5*cRef (robust medians)
%   3) Left edge extends as far as possible while staying fully inside the
%      COMMON valid mask (intersection over datasets), for that vertical band.
%
% Coordinates:
%   X = d.v1 (mm): varies down rows  -> x
%   Z = d.v3 (mm): varies across cols-> z
%
% Output struct rect contains xL,xR,zB,zT and index bounds.

if nargin < 2, nErode = 3; end

files = dir(fullfile(dataFolder,'*.mat'));
assert(~isempty(files), 'No .mat files found in %s', dataFolder);

% --- Reference grid
S0 = load(fullfile(files(1).folder, files(1).name));
d0 = S0.datasa;
X0 = d0.v1;          % x (mm)
Z0 = d0.v3;          % z (mm)

xline = X0(:,1);     % Nrow x 1
zline = Z0(1,:);     % 1 x Ncol

common = true(size(X0));

xR_all  = zeros(numel(files),1);
zLE_all = zeros(numel(files),1);
c_all   = zeros(numel(files),1);

% --- Build common valid mask + collect LE/chord info
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name));
    d = S.datasa;

    valid = (d.v44 == 1);
    valid = keepLargestCC(valid);
    if nErode > 0
        valid = bwmorph(valid,'erode',nErode);
    end
    common = common & valid;

    % --- Leading edge and chord from airfoilrot (mm)
    % Assumed ordering: airfoilrot = [x; z] with LE at min(x).
    % If your airfoilrot is [z; x], swap the two lines below.
    x_air = d.airfoilrot(1,:);
    z_air = d.airfoilrot(2,:);

    [~, iLE] = min(x_air);
    xLE = x_air(iLE);
    zLE = z_air(iLE);

    % If you have chord in the file, use it here. Otherwise keep fixed:
    c = 1;

    c_all(k)   = c;
    zLE_all(k) = zLE;
    xR_all(k)  = xLE - 0.5*c;
end

% --- Conservative right edge that works for all datasets
xR_req = min(xR_all);

% --- Robust vertical placement (lab-fixed)
cRef   = median(c_all);
zLEref = median(zLE_all);

zB = zLEref - 0.5*cRef;
zT = zLEref + 0.5*cRef;

% --- Convert vertical bounds to column indices
jB = find(zline >= zB, 1, 'first');
jT = find(zline <= zT, 1, 'last');
if isempty(jB) || isempty(jT) || jB >= jT
    error('Vertical bounds out of grid. Adjust chord scaling or check Z range.');
end
cols = jB:jT;

% --- Determine x direction
isInc = all(diff(xline) > 0);

% --- Row validity for the vertical band (must be fully valid across cols)
rowGood = all(common(:, cols), 2);

% --- Find iR from requested xR, then snap to a "good" row
if isInc
    % increasing x with row index
    iR0 = find(xline <= xR_req, 1, 'last');
    if isempty(iR0)
        error('xR is outside grid. Check x direction / sign.');
    end

    % choose last good row up to iR0
    iR = find(rowGood(1:iR0), 1, 'last');
    if isempty(iR)
        error(['No valid row found up to requested xR for the chosen vertical band. ', ...
               'Try smaller height (c) or less erosion (nErode).']);
    end

    % left edge = start of contiguous TRUE-run ending at iR
    lastBad = find(~rowGood(1:iR), 1, 'last');
    if isempty(lastBad)
        iL = 1;
    else
        iL = lastBad + 1;
    end

else
    % decreasing x with row index
    iR0 = find(xline <= xR_req, 1, 'first');
    if isempty(iR0)
        error('xR is outside grid. Check x direction / sign.');
    end

    % choose first good row from iR0 downward
    rel = find(rowGood(iR0:end), 1, 'first');
    if isempty(rel)
        error(['No valid row found from requested xR for the chosen vertical band. ', ...
               'Try smaller height (c) or less erosion (nErode).']);
    end
    iR = iR0 + rel - 1;

    % left edge = end of contiguous TRUE-run ending at iR (toward bottom)
    firstBad = find(~rowGood(iR:end), 1, 'first');
    if isempty(firstBad)
        iL = size(common,1);
    else
        iL = iR + firstBad - 2;
    end
end

% --- Rectangle row bounds
i1 = min(iL, iR);
i2 = max(iL, iR);

% --- Snap physical coords to actual used rows
xL = xline(iL);
xR = xline(iR);

% --- Fill ratio check
rectArea = false(size(common));
rectArea(i1:i2, jB:jT) = true;
fillRatio = nnz(rectArea & common) / nnz(rectArea);

% --- Output
rect = struct();
rect.xL = xL; rect.xR = xR;
rect.zB = zB; rect.zT = zT;
rect.i1 = i1; rect.i2 = i2;
rect.j1 = jB; rect.j2 = jT;
rect.cRef = cRef;
rect.zLEref = zLEref;
rect.fillRatio = fillRatio;
rect.nErode = nErode;
rect.xR_req = xR_req;

end

function maskOut = keepLargestCC(maskIn)
cc = bwconncomp(maskIn, 4);
if cc.NumObjects == 0
    maskOut = false(size(maskIn));
    return
end
np = cellfun(@numel, cc.PixelIdxList);
[~, iBig] = max(np);
maskOut = false(size(maskIn));
maskOut(cc.PixelIdxList{iBig}) = true;
end