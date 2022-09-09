function x = newton_raphson_binary_search(f, x_below, x_upper, tol)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
% x  = (x_upper + x_below)/2;
% f_at_x = f(x);
% last_x = x - 1;
% 
% while(abs(f_at_x) > tol && abs(x - last_x) > tol*abs(x))
%     df_at_x = df(x);
%     dx      = f_at_x / df_at_x;
% 
%     last_x = x;
% 
%     if(sign(f_at_x) * sign(df_at_x) < 0)
%         x_below = last_x;
%     else
%         x_upper = last_x;
%     end
% 
%     if(last_x - dx > x_upper || last_x - dx < x_below)
%         % binary search
%         x = (x_below + x_upper)/2;
%     elseif(abs(dx) < (x_upper - x_below)/2)
%         % too little distance
%         x = (x_below + x_upper)/2;
%     else
%         % newton raphson method
%         x = last_x - dx;
%     end
% 
%     f_at_x  = f(x);
% end

f_upper = f(x_upper);
f_below = f(x_below);

fsig = sign(f_upper - f_below);
xsig = sign(x_upper - x_below);

if(sign(f_upper * f_below) ~=-1)
    if(abs(f_upper) < abs(f_below))
        x = x_upper;
        return;
    else
        x = x_below;
        return;
    end
end

diff      = 0;
diff_last = 1;

while(diff < diff_last)
    diff_last = abs(x_upper - x_below);

    x_mid   = x_below + (x_upper - x_below)/2;
    f_mid   = f(x_mid);
    
    x1 = (x_mid*f_below - x_below*f_mid)/(f_below-f_mid);
    fx1 = f(x1);
    x2 = (x_mid*f_upper - x_upper*f_mid)/(f_upper-f_mid);
    fx2 = f(x2);
    
    if(sign(f_mid) == fsig)
        x_upper = x_mid;
        f_upper = f_mid;
    else
        x_below = x_mid;
        f_below = f_mid;
    end

    if(sign(fx1) == fsig && xsig*x1 < xsig*x_upper && xsig*x1 > xsig*x_below)
        x_upper = x1;
        f_upper = fx1;
    end

    if(sign(fx2) == fsig && xsig*x2 < xsig*x_upper && xsig*x2 > xsig*x_below)
        x_upper = x2;
        f_upper = fx2;
    end

    if(sign(fx1) ~= fsig && xsig*x1 > xsig*x_below && xsig*x1 < xsig*x_upper)
        x_below = x1;
        f_below = fx1;
    end

    if(sign(fx2) ~= fsig && xsig*x2 > xsig*x_below && xsig*x2 < xsig*x_upper)
        x_below = x2;
        f_below = fx2;
    end

    diff = abs(x_upper - x_below);
end

x = (x_below + x_upper)/2;

end