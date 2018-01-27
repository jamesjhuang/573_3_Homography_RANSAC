function [bestH, bestPairL, bestPairR, residual] = ...
    ransac_est(pairLeft, pairRight, iterations, model, numRand, norm)

    iter = iterations;
    inliers = zeros(iter, 1);
    for i = 1:iter
        % Selecting random pairs, here row = x and col = y
        randPairs = randperm(length(pairLeft), numRand);
        xL = pairLeft(randPairs, 1);  % x2
        yL = pairLeft(randPairs, 2);  % y2
        xR = pairRight(randPairs, 1); % x1
        yR = pairRight(randPairs, 2); % y1

        if norm == 1
            matches = [xL yL xR yR];
        else
            matches = [xR yR xL yL]; 
        end
        
        h = model(matches, norm);

        if norm == 0
            % Using obtained homography matrix to check for number of inliers
            for j = 1:length(pairRight)
                % Equation is of the form xH=x', in this case xH=b
                % So b represents the new coordinates
                b = h * [pairRight(j,1); pairRight(j,2); 1]; 
                estimatedCoord(j, 1:2) = b(1:2) / b(3);
                residual(j) = sqrt((estimatedCoord(j,1) - pairLeft(j, 1))^2 ...
                    + (estimatedCoord(j,2) - pairLeft(j, 2))^2);
            end


        else
            N = length(pairRight);
            matches = [pairLeft pairRight];
            L = (h * [matches(:,1:2) ones(N,1)]')'; % transform points from the first
                                    % image to get epipolar lines in the second image
            % find points on epipolar lines L closest to matches(:,3:4)
            L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
            residual = sum(L .* [matches(:,3:4) ones(N,1)], 2);
            estimatedCoord = pairLeft;
        end
        
        inliers(i) = numel(find(residual <= 10));
        inNout(i) = inliers(i) / length(pairRight);

        if inNout(i) > 0.4
            inlierRow = find(residual <= 10);
            x1 = pairRight(inlierRow, 1);
            y1 = pairRight(inlierRow, 2);
            x2 = pairLeft(inlierRow, 1);
            y2 = pairLeft(inlierRow, 2);
            matches = [x1 y1 x2 y2];
            h = model(matches, 1);
            if inliers(i) == max(inliers)
                bestH = h;
                bestPairL = [x2 y2]; % [x y]
                bestPairR = [x1 y1];      % [x y]
            end
        end

        
    end
    max(inliers)
end