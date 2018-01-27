% CSE573
% Written by: James J. Huang


function homography = Homography(matches, norm)
    useless = norm;
    x1 = matches(:,1); 
    y1 = matches(:,2);
    x2 = matches(:,3);
    y2 = matches(:,4);
    
    % Forming matrix to solve nonhomogeneous system of the form xH=x'
    j = 1;
     for i = 1:length(x1)
        x(j:j+1, :) = [x1(i) y1(i) 1 0 0 0 -x1(i)*x2(i) -y1(i)*x2(i);
                       0 0 0 x1(i) y1(i) 1 -x1(i)*y2(i) -y1(i)*y2(i)];
        j = j+2;
     end
     x_p(1:2:2*length(x1)-1, 1) = x2(:);
     x_p(2:2:2*length(x1), 1) = y2(:);
     
    % Solving for homography matrix
    h = x \ x_p;
    % Above nonhomogeneous system assumes h(3,3) = 1
    h(9) = 1;
    % Reshaping column of h into 3x3 matrix
    homography = reshape(h, [3 3])';
%     [~, ~, V] = svd(x);
%     homography = reshape(V(:,end), [3 3]);
%     homography = homography ./ homography(3,3);
end