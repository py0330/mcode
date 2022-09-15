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

Tmax   = zeros(1,n);
Tmin   = zeros(1,n);


if(T == inf)
    T=1000;
end

% 因为输入的 T 对本次规划肯定是成功的，因此如果是最后一个周期，那么就直接成功
if(i == m)
    success = 1;
    return;
else
   for j=1:n
        [v0(j),~,~,~]      = s_make_s_curve( ...
            p0(j), ...
            v0(j), ...
            pos(i,j), ...
            vel(i,j), ...
            acc(i,j), ...
            jerk(i,j), ...
            T);
        p0(j) = pos(i,j);
        
        [Tmax(j), Tmin(j)] = s_s_curve_Tmax_Tmin( ...
            p0(j), ...
            v0(j), ...
            pos(i+1,j), ...
            max_vb(i+1,j), ...
            vel(i+1,j), ...
            acc(i+1,j), ...
            jerk(i+1,j));
    end

    Tmax_all = min(Tmax(:));
    Tmin_all = max(Tmin(:));

    if(Tmax_all == inf)
%         fprintf('test slow success:%d \n T:%f \n Tmax_ALL:%f \n Tmin_all:%f',i,T,Tmax_all,Tmin_all)
%         Tmax
%         Tmin
%         v0
        success = 1;
        return;
    end

    if(Tmax_all == -1 || Tmax_all < Tmin_all)
%         fprintf('test slow failed:%d \n T:%f \n Tmax_ALL:%f \n Tmin_all:%f',i,T,Tmax_all,Tmin_all)
%         Tmax
%         Tmin
%         v0
        success = 0;
        return;
    end
    
    success = s_test_curve_slow(i + 1, p0, v0, pos, max_vb, vel, acc, jerk, Tmax_all);
    return;

end



end


