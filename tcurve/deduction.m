%%
syms pb pe vb ve vmax amax pt v a T Ta Tb
%%

eq = Ta*(vb+v)/2 + Tb*(ve+v)/2 + (T-Ta-Tb)*v
collect(eq,v)

Ta = (v-vb)/amax
Tb = (v-ve)/amax
eq2 = Ta*(vb+v)/2 + Tb*(ve+v)/2
eq21 = Ta*(vb+v)/2 + Tb*(ve+v)/2 + (T-Ta-Tb)*v
collect(eq2,v)

Ta = (-v+vb)/amax
Tb = (-v+ve)/amax
eq3 = Ta*(vb+v)/2 + Tb*(ve+v)/2
eq31 = Ta*(vb+v)/2 + Tb*(ve+v)/2 + (T-Ta-Tb)*v
collect(eq3,v)

%%

% vb = 1.5;
% ve = 1.5;
% pb = 0;
% pe = 0.1;
% vmax = 3;
% amax = 1.5;

vb = 1.2412986333062634;
ve = 1.2838800173112299;
pb = 0.052814220624860846;
pe = 0.055368645958708947;
vmax = 3.1400000000000001;
amax = 31.399999999999999;

%%
[T1,T2,T3] = s_tcurve_T_range(pb,pe,vb,ve,vmax,amax);
T=0.295;
[Ta,Tb,v,a,mode] = s_tcurve_param(pb, pe, vb, ve, vmax, amax, T);

dt = 0.001;
t = 0:dt:T;




p = zeros(size(t));
for i =1:length(t)
    p(i) = s_tcurve_value(pb, pe, vb, ve, T, Ta, Tb, mode, v, a, t(i));
end

%% plot
subplot(1,3,1);
plot(t,p)
subplot(1,3,2)
plot(diff(p)/dt)
subplot(1,3,3)
plot(diff(diff(p))/dt/dt)


%% 
n=1e8;
for i=1:n
    i
    pb = (rand - 0.5)/(rand+0.05);
    pe = (rand - 0.5)/(rand+0.05);
    vmax = (rand + 0.1)/(rand+0.05);
    amax = (rand + 0.1)/(rand+0.05);
    vb = (rand - 0.5)/(rand+0.05);
    vb = max(vb,-vmax);
    vb = min(vb,vmax);
    ve = (rand - 0.5)/(rand+0.05);
    ve = max(ve,-vmax);
    ve = min(ve,vmax);
    
    [T1,T2,T3] = s_tcurve_T_range(pb,pe,vb,ve,vmax,amax);
    
    for T = [T1:(T2-T1)/10:T2, T3:T3/10:2*T3]
        [Ta,Tb,v,a,mode] = s_tcurve_param(pb, pe, vb, ve, vmax, amax, T);

        t = linspace(0,T,1000);
        dt =t(2)-t(1);
        
        p = zeros(size(t));
        for i =1:length(t)
            p(i) = s_tcurve_value(pb, pe, vb, ve, T, Ta, Tb, mode, v, a, t(i));
        end
        
        if(max(diff(p))/dt > vmax + 1e-12 / dt || min(diff(p))/dt < -vmax-1e-12 / dt)
            error("v error")
        end

        if(max(diff(diff(p)))/dt/dt > amax + max(2e-14/dt/dt,1e-8) || min(diff(diff(p)))/dt/dt < -amax -max(2e-14/dt/dt,1e-8))
            error("a error")
        end

        
        if(abs(p(1) - pb)>1e-10)
            error("pb error")
        end

        if(abs(p(end) - pe)>1e-10)
            error("pe error")
        end

    end
    
    



end















