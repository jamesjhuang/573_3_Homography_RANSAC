% CSE593 Homework 3
% Written by: James Huang
% Date: 10-26-2017
% Main driver script


%% 1. Stitching 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~ Part I. Stitching ~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Loading left and right images
imgL = imread('../data/part1/uttower/left.jpg');
imgL = im2double(imgL);
imgR = imread('../data/part1/uttower/right.jpg');
imgR = im2double(imgR);

[H, W, Z] = size(imgL);

% Converting to grayscale then to double
imgLG = rgb2gray(imgL);
imgLG = im2double(imgLG);
imgRG = rgb2gray(imgR);
imgRG = im2double(imgRG);


% Detecting features using harris.m
radius = 2;
sigma = sqrt(2);
thresh = 0.1;
disp = 1;
neighborhood = 31;

% Using harris to detector corners
[cimL, rL, cL] = harris(imgLG, sigma, thresh, radius, disp);
[cimR, rR, cR] = harris(imgRG, sigma, thresh, radius, disp);

% Padding image w/ replicate pixels
imgLGP = padarray(imgLG, [neighborhood neighborhood], 'replicate');
imgRGP = padarray(imgRG, [neighborhood neighborhood], 'replicate');

% Shifting the row and column index due to padding
rL = rL + neighborhood;
cL = cL + neighborhood;
rR = rR + neighborhood;
cR = cR + neighborhood;

% Extract local neighborhoods and flattening them
neigh = (neighborhood - 1) / 2;
descriptorL = zeros(length(rL), neighborhood^2);
descriptorR = zeros(length(rR), neighborhood^2);

for i = 1:length(rL)
    box = imgLGP(rL(i) - neigh:rL(i) + neigh, ...
                cL(i) - neigh:cL(i) + neigh);
    % Normalizing for 0 mean and unit standard deviation for each
    % descriptor independently
    descriptorL(i, :) = zscore(box(:))';
end

for i = 1:length(rR)
    box = imgRGP(rR(i) - neigh: rR(i) + neigh, ...
        cR(i) - neigh: cR(i) + neigh);
    descriptorR(i, :) = zscore(box(:))';
end

% Distance between descriptors of 2 images
desDist = dist2(descriptorL, descriptorR);

% Selecting putative matches using smallest x pairs
numPairs = 100;
[pairR, pairC] = putativeMatches(desDist, numPairs);

% Shifting the row and column index due to padding
rL = rL - neighborhood;
cL = cL - neighborhood;
rR = rR - neighborhood;
cR = cR - neighborhood;

% Extracting coordinates of matches on image
% pairR and pairC obtained from above represents the descriptor, not the
% actual coordinates. So, we use these descriptors to extract the
% coordinates. According to dist2, pairR will represent descriptors from
% the L image, and pairC represents descripts from the R image.
pairLeft = [cL(pairR) rL(pairR)];
pairRight = [cR(pairC) rR(pairC)];

% Plotting pairs that passed
plotPairs(imgLG, imgRG, pairLeft, pairRight);

% RANSAC for estimating homography matrix
[bestH, bestPairL, bestPairR, residual] = ...
    ransac_est(pairLeft, pairRight, 5000, @Homography, 4, 0);

% Refining pairs by removing those out of bounds
% badR = find((bestPairL(:,1) >= W));
% bestPairL(badR,:) = [];
% bestPairR(badR,:) = [];
% badC = find((bestPairL(:,2) >= H));
% bestPairL(badC,:) = [];
% bestPairR(badC,:) = [];

% Plotting pairs again
% plotPairs(imgLG, imgRG, bestPairL, bestPairR);

T = maketform('projective', inv(bestH'));

[imgLT, xLims, yLims] = imtransform(imgL, T);
yLims = abs(ceil(yLims));
xLims = abs(ceil(xLims));
pano = zeros(H + yLims(1), W + xLims(1), 3);
[Ht, Wt, Zt] = size(pano);
figure;
for i = 1:3
    stitch = zeros(H + yLims(1), W + xLims(1), 3);
    stitch(Ht-H+1:Ht, Wt-W+1:Wt, 2) = imgR(:,:,i);
    stitch(1:size(imgLT, 1), 1:size(imgLT, 2), 1) = imgLT(:,:,i);
    stitch(stitch(:,:,:) == 0) = NaN;
    stitch(:,:,3) = mean(stitch,3, 'omitnan');
    stitch(isnan(stitch)) = 0;
    pano(:,:,i) = stitch(:,:,3);
end
panoRGB = cat(3, pano(:,:,1), pano(:,:,2), pano(:,:,3));
imshow(panoRGB);


%% 2. Fundamental Matrix and Triangulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~ Part II. Fun Matrix and Triangulation ~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% load images and match files for the first example
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

I1 = imread('../data/part2/library1.jpg');
I2 = imread('../data/part2/library2.jpg');
matches = load('../data/part2/house_matches.txt'); 
% this is a N x 4 file where the first two numbers of each row
% are coordinates of corners in the first image and the last two
% are coordinates of corresponding corners in the second image: 
% matches(i,1:2) is a point in the first image [x y] [c r]
% matches(i,3:4) is a corresponding point in the second image

N = size(matches,1);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% display two images side-by-side with matches
% this code is to help you visualize the matches, you don't need
% to use it to produce the results for the assignment
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% display second image with epipolar lines reprojected 
% from the first image
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% first, fit fundamental matrix to the matches
norm = 1;
F = fit_fundamental(matches, norm); % this is a function that you should write
L = (F * [matches(:, 1:2) ones(N,1)]')'; % transform points from the first
                        % image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)], 2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

mean(pt_line_dist)

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
figure;
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

F = fit_fundamental(matches, norm); % this is a function that you should write
L = (F * [matches(:, 1:2) ones(N,1)]')'; % transform points from the first
                        % image to get epipolar lines in the second image

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,1:2) ones(N,1)], 2);
closest_pt = matches(:,1:2) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

mean(pt_line_dist)

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
figure;
imshow(I1); hold on;
plot(matches(:,1), matches(:,2), '+r');
line([matches(:,1) closest_pt(:,1)]', [matches(:,2) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~ RANSAC Implementation ~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

I1G = im2double(rgb2gray(I1));
I2G = im2double(rgb2gray(I2));

% Detecting features using harris.m
radius = 2;
sigma = sqrt(2);
thresh = 0.005;
disp = 1;
neighborhood = 3;

% Using harris to detector corners
[cimL, rL, cL] = harris(I1G, sigma, thresh, radius, disp);
[cimR, rR, cR] = harris(I2G, sigma, thresh, radius, disp);
%%
% Padding image w/ replicate pixels
I1G = padarray(I1G, [neighborhood neighborhood], 'replicate');
I2G = padarray(I2G, [neighborhood neighborhood], 'replicate');

% Shifting the row and column index due to padding
rL = rL + neighborhood;
orhood;cL = cL + neighborhood;
rR = rR + neighborhood;
cR = cR + neighb

% Extract local neighborhoods and flattening them
neigh = (neighborhood - 1) / 2;
descriptorL = zeros(length(rL), neighborhood^2);
descriptorR = zeros(length(rR), neighborhood^2);

for i = 1:length(rL)
    box = I1G(rL(i) - neigh:rL(i) + neigh, ...
                cL(i) - neigh:cL(i) + neigh);
    % Normalizing for 0 mean and unit standard deviation for each
    % descriptor independently
    descriptorL(i, :) = (box(:))';
end

for i = 1:length(rR)
    box = I2G(rR(i) - neigh: rR(i) + neigh, ...
        cR(i) - neigh: cR(i) + neigh);
    descriptorR(i, :) = (box(:))';
end

% Distance between descriptors of 2 images
desDist = dist2(descriptorL, descriptorR);

% Selecting putative matches using smallest x pairs
numPairs = 75;
[pairR, pairC] = putativeMatches(desDist, numPairs);

% Shifting the row and column index due to padding
rL = rL - neighborhood;
cL = cL - neighborhood;
rR = rR - neighborhood;
cR = cR - neighborhood;

% Extracting coordinates of matches on image
% pairR and pairC obtained from above represents the descriptor, not the
% actual coordinates. So, we use these descriptors to extract the
% coordinates. According to dist2, pairR will represent descriptors from
% the L image, and pairC represents descripts from the R image.
pairLeft = [cL(pairR) rL(pairR)];
pairRight = [cR(pairC) rR(pairC)];

% Plotting pairs that passed
plotPairs(I1, I2, pairLeft, pairRight);

%%
norm = 1;
[bestF, bestPairL, bestPairR, residual] = ...
    ransac_est(pairLeft, pairRight, 5000, @fit_fundamental, 8, norm);


L = (bestF * [matches(:,1:2) ones(N,1)]')'; % transform points from the first
                        % image to get epipolar lines in the second image

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)], 2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

mean(pt_line_dist)

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
figure;
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% find points on epipolar lines L closest to matches(:,3:4)
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,1:2) ones(N,1)], 2);
closest_pt = matches(:,1:2) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

mean(pt_line_dist)

% find endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

% display points and segments of corresponding epipolar lines
figure;
imshow(I1); hold on;
plot(matches(:,1), matches(:,2), '+r');
line([matches(:,1) closest_pt(:,1)]', [matches(:,2) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');



%%
P1 = load('../data/part2/house1_camera.txt');
P2 = load('../data/part2/house2_camera.txt');

% Finding camera centers using SVD, output is in homogeneous 
% coordinates so divide through by last element to recover the 
% cartesian coordinates
[~, ~, V] = svd(P1);
center1 = V(:,end);
center1 = center1 ./ center1(end);

[~, ~, V] = svd(P2);
center2 = V(:,end);
center2 = center2 ./ center2(end);


P11 = load('../data/part2/library1_camera.txt');
P22 = load('../data/part2/library2_camera.txt');

[~, ~, V] = svd(P11);
center11 = V(:,end);
center11 = center11 ./ center11(end);

[~, ~, V] = svd(P22);
center22 = V(:,end);
center22 = center22 ./ center22(end);

matches = load('../data/part2/house_matches.txt'); 

% Converting ground truth matches to homogeneous by appending 1
coord1 = matches(:,1:2);
coord1(:,3) = 1;
coord2 = matches(:,3:4);
coord2(:,3) = 1;

for i = 1:length(matches)
    % Using the linear approach where [x_1x]P1X = 0 and
    % [x_2x]P2X = 0 to solve for X, which represents the coordinates
    % of the match in 3D coordinates
    x_1x = [0 -coord1(3) coord1(2); coord1(3) 0 -coord1(1);
            -coord1(2) coord1(1) 0];
    x_2x = [0 -coord2(3) coord2(2); coord2(3) 0 -coord2(1);
            -coord2(2) coord2(1) 0];
    % homogeneous1/2 are 3x4 matrices, representing homogeneous coordinates
    homogeneous1 = x_1x * P1;
    homogeneous2 = x_2x * P2;

    A(1,:) = [coord1(i,2) * P1(3,:) - P1(2,:)];
    A(2,:) = [P1(1,:) - coord1(i,1) * P1(3,:)];
    A(3,:) = [coord2(i,2) * P2(3,:) - P2(2,:)];
    A(4,:) = [P2(1,:) - coord2(i,1) * P2(3,:)];

    [U, S, V] = svd(A);
    coord3(i, :) = V(:,end) ./ V(end, end);
end

figure;
plot3(coord3(:,1), coord3(:,2), coord3(:,3), 'o');
axis equal; grid on; hold on;
plot3(center1(1), center1(2), center1(3), 'ro');
plot3(center2(1), center2(2), center2(3), 'go');
legend('X', 'Cam1', 'Cam2');
title('House');
