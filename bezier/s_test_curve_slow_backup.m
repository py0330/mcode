function success = s_test_curve_slow(i, p0, v0, pos, max_vb, vel, acc, jerk, T)
% 尝试将速度降下来
%
% pa     : current pos
% va     : current vel
% pos    : target positions
% max_vb : max vel at target positions
% vel    : max vel  during period
% acc    : max acc  during period
% jerk   : max jerk during period

[m,n] = size(pos);

Tmax   = zeros(m,n);
Tmin   = zeros(m,n);

for j = 1:n
    [Tmax(i,j), Tmin(i,j)] = s_s_curve_Tmax_Tmin( ...
        p0(j), ...
        v0(j), ...
        pos(i,j), ...
        max_vb(i,j), ...
        vel(i,j), ...
        acc(i,j), ...
        jerk(i,j));

    if(Tmax(i,j) == -1)
        success = 0;
        return;
    end
end

Tmax_all = min(Tmax(i,:));
Tmin_all = max(Tmin(i,:));

if(Tmax_all == inf)% 可以停下来
    success = 1;
    return;
elseif(Tmin_all > Tmax_all)%无法规划
    success = 0;
    return;
elseif(i == m)
    success = 1;
    return;
else
    for j=1:n
        [v0(j),~,~,~]=s_make_s_curve( ...
            p0(j), ...
            v0(j), ...
            pos(i+1,j), ...
            vel(i+1,j), ...
            acc(i+1,j), ...
            jerk(i+1,j), ...
            Tmax_all);
    end

    p0 = pos(i+1,:);
    success = s_test_curve_slow(i + 1, p0, v0, pos, max_vb, vel, acc, jerk);
    return;

end



end


