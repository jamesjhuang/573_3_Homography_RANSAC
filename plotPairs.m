function plotPairs(imgL, imgR, leftCoord, rightCoord)
    figure;
    imshowpair(imgL,imgR,'montage');
    hold on;
    rightCoord(:,1) = rightCoord(:,1) + size(imgL, 2);
    for i = 1:length(rightCoord)
        plot([rightCoord(i,1) leftCoord(i,1)], ...
            [rightCoord(i,2) leftCoord(i,2)], 'r', 'linewidth', 1);
        hold on;
        plot([rightCoord(i,1) leftCoord(i,1)], ...
            [rightCoord(i,2) leftCoord(i,2)], 'yo', 'linewidth', 1);
        hold on;
    end
end