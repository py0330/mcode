%%
clear;
pa = 0;
pb = 0.8;
vb = 1.5;
va_upper = 1.6;
va_below = 0;
vc_max = 2;
vb_max = 0.5;
a = 5;
j=10;
T=1;

va_range = va_below:0.001:va_upper;

sz = size(va_range);
vb_uppers = zeros(sz);
vb_belows = zeros(sz);
for(i=1:sz(2))
    if(va_range(i)>=1.0)
        112
%         message('debug')
    end
    va = va_range(i);
    
    vb_uppers(i) = s_scurve_cpt_vb_upper(pa,va,pb,vc_max,vb_max,a,j,T);
    vb_belows(i) = s_scurve_cpt_vb_below(pa,va,pb,vc_max,vb_max,a,j,T);
    
    

end

plot(vb_uppers);

hold on
plot(vb_belows);
hold on
[vb_range_upper, vb_range_below] = s_scurve_cpt_vb_range(pa,pb, 1.6, 0,vc_max,vb_max,a,j,T);

line([0,4000], [vb_range_upper,vb_range_upper],'Color','k','LineStyle','--');
hold on
line([0,4000], [vb_range_below,vb_range_below],'Color','r','LineStyle','--');


