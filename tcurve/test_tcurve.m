function test_tcurve(n)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明

for i=1:n
    vb = (rand - 0.5)/(rand+0.05);
    ve = (rand - 0.5)/(rand+0.05);
    pb = (rand - 0.5)/(rand+0.05);
    pe = (rand - 0.5)/(rand+0.05);
    vmax = (rand + 0.1)/(rand+0.05);
    amax = (rand + 0.1)/(rand+0.05);
    
    [T1,T2,T3] = s_tcurve_T_range(pb,pe,vb,ve,vmax,amax);
    
    for T = [T1:(T2-T1)/10:T2, T3:T3/10:2*T3]
        [Ta,Tb,v,a,mode] = s_tcurve_param(pb, pe, vb, ve, vmax, amax, T);

        


%         dt = T/1000;
        t = linspace(0,T,1000);
        dt =t(2)-t(1);
        
        p = zeros(size(t));
        for i =1:length(t)
            p(i) = s_tcurve_value(pb, pe, vb, ve, T, Ta, Tb, mode, v, a, t(i));
        end
        
        if(max(diff(p))/dt > vmax+1e-10 || min(diff(p))/dt < -vmax-1e-10)
            error("v error")
        end

        if(max(diff(diff(p)))/dt/dt > amax+1e-10 || max(diff(diff(p)))/dt < -amax-1e-10)
            error("a error")
        end

        
        if(abs(p(1) - pb)<1e-10)
            error("pb error")
        end

        if(abs(p(end) - pe)<1e-10)
            error("pe error")
        end

    end
    
    



end



end