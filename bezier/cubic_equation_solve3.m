function x = cubic_equation_solve3(a,b,c,d)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

% %%%%%%%%%  solve method 1 %%%%%%%%%%
% p = (3*a*c - b^2) / (3*a^2);
% q = (2*b^3-9*a*b*c + 27*a^2*d) / (27*a^3);
% 
% C0 = (-q/2+sqrt(q^2/4+p^3/27))^(1/3);
% C1 = (-q/2+sqrt(q^2/4+p^3/27))^(1/3) * (-0.5 - sqrt(3)/2*1i);
% C2 = (-q/2+sqrt(q^2/4+p^3/27))^(1/3) * (-0.5 + sqrt(3)/2*1i);
% 
% C = [C0;C1;C2];
% 
% y = C - p ./ (3 * C);
% 
% x = y - b/(3*a);
% 
% for k =1:length(x)
%     if(abs(imag(x(k))) < 1e-10)
%         x(k) = real(x(k));
%     end
% end

% %%%%%%%%%  solve method 2 %%%%%%%%%%
% % see https://stackoverflow.com/questions/13328676/c-solving-cubic-equations
% %
% % 以下当b\c很小时，无法求得准确结果，例如 b c接近1e-6，此时disc中的q接近
% b = b/a;
% c = c/a;
% d = d/a;
% 
% q = (3.0*c - (b*b))/9.0;
% r = (-(27.0*d) + b*(9.0*c - 2.0*(b*b)))/54.0;
% 
% x = zeros(1,3);
% 
% % disc = q*q*q + r*r;
% disc = (b^3*d)/27 - (b^2*c^2)/108 - (b*c*d)/6 + c^3/27 + d^2/4;
% term1 = (b/3.0);
% if(disc > 0)
%     s = r + sqrt(disc);
%     s = sign(s)*abs(s)^(1/3);
%     t = r - sqrt(disc);
%     t = sign(t)*abs(t)^(1/3);
%     
%     x(1) = -term1 + s + t;
%     term1 = term1 + (s + t)/2.0;
%     x2_real = -term1;
%     x3_real = -term1;
%     term1 = sqrt(3.0)*(-t + s)/2;
%     x2_imag = term1;
%     x3_imag = -term1;
%     x(2) = complex(x2_real, x2_imag);
%     x(3) = complex(x3_real, x3_imag);
% elseif(disc == 0)
%     r13 = sign(r)*abs(r)^(1/3);
%     x(1) = -term1 + 2.0*r13;
%     x(2) = -(r13 + term1);
%     x(3) = x(2);
%     x = sort(x);
% else
%     q=-q;
%     dum1 = q*q*q;
%     dum1 = acos(r/sqrt(dum1));
%     r13 = 2.0*sqrt(q);
%     x(1) = -term1 + r13*cos(dum1/3.0);
%     x(2) = -term1 + r13*cos((dum1 + 2.0*pi)/3.0);
%     x(3) = -term1 + r13*cos((dum1 + 4.0*pi)/3.0);
% 
%     x = sort(x);
% end

%%%%%%%%%  solve method 3 %%%%%%%%%%
% 对于n次方程：
%
% f(x) = x^n + kn-1 * x^n-1 + kn-2 * x^n-2 + ... + k1 * x + k0
%
% 若 x = 1 + |k0| + |k1| ... |kn-1|
% f(x) = (1 + |k0| + |k1| ... |kn-1|)*x^n-1 + kn-1 * x^n-1 + kn-2 * x^n-2 + ... + k1 * x + k0
%
% 展开并进行结合，有：
% f(x) = [x^n-1 * (1 + |k0|) + k0 ] + [x*(x^n-2 * |k1| + k1)] +...+ [x^n-1*(|kn-1| + kn-1)]
%      > 0
% 
% 同理，若x = -1 - |k0| - |k1| ... |kn-1|
% 则必有 f(x) < 0
%
% --------------------------------------------------
% 对于3次方程 x^3 + b*x^2 + c*x + d = 0  
% diff(x^3 + b*x^2 + c*x + d)
% = 3*x^2 + 2*b*x + c
%
% 两个极值点为：
% xp1 = (-b - sqrt(b*b - 3*c)) / 3
% xp2 = (-b + sqrt(b*b - 3*c)) / 3

b = b/a;
c = c/a;
d = d/a;

f  = @(x)local_f(x,b,c,d);

x_lhs = -1 - abs(b) - abs(c) - abs(d);
x_rhs = -x_lhs;

% CASE 0: 没有极点，b*b - 3*c < 0
if(b*b - 3*c < 0)
    x = zeros(1,1);
    x_below = x_lhs;
    x_upper = x_rhs;
    x(1) = newtor_raphson_binary_search(f, x_upper,x_below,10*eps);
    return;
end

xp1 = (-b - sqrt(b*b - 3*c)) / 3;
xp2 = (-b + sqrt(b*b - 3*c)) / 3;

f1    = f(xp1);
f2    = f(xp2);



% CASE 1: xp1 == xp2 且 f(xp1) == f(xp2) == 0
% 此时有3重根：
% x1 = x2 = x3 = xp1 == xp2
if(abs(f1) < eps && abs(f2) < eps)
    x = zeros(1,3);
    x(1) = (xp1 + xp2)/2;
    x(2) = x(1);
    x(3) = x(1);
    return;
end
% CASE 2: f(xp1) == 0
% 此时有个2重根和1重根：
%  (-1-|b|-|c|-|d|) < x1 = x2 < xp1
%              xp2  < x3      < 1+|b|+|c|+|d|
if(abs(f1) < eps)
    x = zeros(1,3);
    x_upper = xp1;
    x_below = x_lhs;
    x(1) = newtor_raphson_binary_search(f, x_upper,x_below,10*eps);
    x(2) = x(1);
    x_upper = x_rhs;
    x_below = xp2;
    x(3) = newtor_raphson_binary_search(f, x_upper,x_below,10*eps);
    return;
end
% CASE 3: f(xp2) == 0
% 此时有个2重根和1重根：
%  (-1-|b|-|c|-|d|) < x1      < xp1
%              xp2  < x2 = x3 < 1+|b|+|c|+|d|
if(abs(f2) < eps)
    x = zeros(1,3);
    x_below = x_lhs;
    x_upper = xp1;
    x(1) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    x_below = xp2;
    x_upper = x_rhs;
    x(2) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    x(3) = x(2);
    return;
end
% CASE 4: f(xp1) * f(xp2) < 0
% 此时有3个根，分别位于：
% (-1-|b|-|c|-|d|) < x1 < xp1
%             xp1  < x2 < xp2
%             xp3  < x3 < 1+|b|+|c|+|d|
if(f1 * f2 < 0)
    x = zeros(1,3);
    x_below = x_lhs;
    x_upper = xp1;
    x(1) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    x_below = xp1;
    x_upper = xp2;
    x(2) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    x_below = xp2;
    x_upper = x_rhs;
    x(3) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    return;
end
% CASE 5: f(xp2) > 0 && f(xp1) > 0
% 此时有1个根，位于：
%       -1-|b|-|c|-|d| < x1 < xp1
if(f1>= 0 && f2 >=0)
    x = zeros(1,1);
    x_below = x_lhs;
    x_upper = xp1;
    x(1) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    return;
end
% CASE 6: f(xp2) < 0 && f(xp1) < 0
% 此时有1个根，位于：
%      xp2  < x1 < 1+|b|+|c|+|d|
if(f1 < 0 && f2 < 0)
    x = zeros(1,1);
    x_below = xp2;
    x_upper = x_rhs;
    x(1) = newton_raphson_binary_search(f, x_upper,x_below,10*eps);
    return;
end    

end % end of main func

function f = local_f(x,b,c,d)
    f = x*(x*(x + b) + c) + d;
    return;
end