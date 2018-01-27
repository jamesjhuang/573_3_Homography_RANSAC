% CSE573
% Written by: James J. Huang

function fundamental = fit_fundamental(matches, norm)
    x1 = matches(:, 1);
    y1 = matches(:, 2);
    x2 = matches(:, 3);
    y2 = matches(:, 4);
    n = length(x1);
    if norm == 1
        % Normalized algorithm
        centroid1 = mean([x1 y1]);
        centroid2 = mean([x2 y2]);
        u = [x1 y1 ones(n, 1)]';
        u_p = [x2 y2 ones(n, 1)]';
        sum1 = 0;
        sum2 = 0;
        for i = 1:n
            sum1 = sum1 + sum(([x1(i) y1(i)] - centroid1).^2);
            sum2 = sum2 + sum(([x2(i) y2(i)] - centroid2).^2);
        end
        s1 = sqrt((0.5/n)*sum1);
        s2 = sqrt((0.5/n)*sum2);
        T1 = [inv(s1) 0 -inv(s1)*centroid1(1);...
              0 inv(s1) -inv(s1)*centroid1(2);...
              0 0 1];
        T2 = [inv(s2) 0 -inv(s2)*centroid2(1);...
              0 inv(s2) -inv(s2)*centroid2(2);...
              0 0 1];
        for i = 1:n
            u(:, i) = T1 * u(:, i);
            u_p(:, i) = T2 * u_p(:, i);
        end
        u = u';
        u_p = u_p';
        x1 = u(:, 1); 
        y1 = u(:, 2);
        x2 = u_p(:,1);
        y2 = u_p(:, 2);
        x = zeros( length(matches), 9);
        for i = 1:length(matches)
            x(i, :) = [x2(i) * x1(i) x2(i) * y1(i) x2(i)...   
                    y2(i) * x1(i) y2(i) * y1(i) y2(i)...
                    x1(i) y1(i) 1];
        end
        [~, ~, V] = svd(x);
        fundamental = reshape(V(:,end), [3 3])';
        fundamental = T2' * fundamental * T1;
    else
        % Unnormalized algorithm
        x = zeros( length(matches), 9);
        for i = 1:length(matches)
            x(i, :) = [x2(i) * x1(i) x2(i) * y1(i) x2(i)...   
                    y2(i) * x1(i) y2(i) * y1(i) y2(i)...
                    x1(i) y1(i) 1];
        end
        [~, ~, V] = svd(x);
        F = reshape(V(:,end), [3 3])';
        [U, S, V] = svd(F);
        S(3,3) = 0;
        fundamental = U * S * V';
    end
end
