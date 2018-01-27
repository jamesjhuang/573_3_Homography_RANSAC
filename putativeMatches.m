% CSE573_HW3
% Written by: James J. Huang

function [rP, cP] = putativeMatches(desDist, numPairs)
    rP = zeros(numPairs, 1);
    cP = zeros(numPairs, 1);
    iter = 1;
    while (rP(numPairs) == 0)
    %for i = 1:numPairs
        % Extracting the row and column coords of min distance returned
        % 'find' in this use will return index of element, so we need to 
        % convert that to subscripts
        [r, c] = ind2sub(size(desDist), find(desDist == min(desDist(:))));
        for j = 1:length(r)
            % In case there is more than 1 pair with the same min distance
            NNsorted = sort(desDist(:, c(j)));
            secNN = find(desDist(:, c(j)) == NNsorted(2));
            NNratio = desDist(r(j), c(j)) / desDist(secNN(1), c(j));
            if NNratio > 0.8
                desDist(r(j), c(j)) = 10000;
                desDist(secNN(1), c(j)) = 10000;
                continue
            else
                rP(iter) = r(j);
                cP(iter) = c(j);
                desDist(r(j), c(j)) = 10000;
                iter = iter+1;
            end
        end
    end
end