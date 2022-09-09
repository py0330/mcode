function x = square_equation_solve(a,b,c)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

x = [(-b+sqrt(b^2-4*a*c))/(2 * a);(-b-sqrt(b^2-4*a*c))/(2 * a)];

for k =1:length(x)
    if(abs(imag(x(k))) < 1e-10)
        x(k) = real(x(k));
    end
end

end