function v = s_acc_vend(va,a,j,T)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
    if(a/j > T/2)
        v = va + j*T*T/4;
    else
        v = va + T*a - a^2/j;
    end
end
