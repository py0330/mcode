%% 一、整体规划规划流程
% 1. 问题描述
% 在一段轨迹中，固定参数如下：
% a :最大加速度
% j :最大加加速度
% pt:位移 
%
% 以上三者基本不可调，一段轨迹中只包含有加速段，匀速段、减速段
% 三种情况，因此对于 a(j) 来说，只能为[a, 0, -a](j) 三种情况。
%
% 其他可调的参数包括：
% va:起始速度
% va_upper:起始速度的上界
% vb_below:起始速度的下界
%
% vb:结束速度
% vb_upper:结束速度的上界
% vb_below:结束速度的下界
% vb_max  :最大可能的结束速度
%
% vc:匀速段速度
% vc_upper:匀速段速度的上界
% vc_below:匀速段速度的下界
% vc_max  :最大可能的匀速段速度
% 
% T :轨迹周期
% T_upper :周期的上限
% T_below :周期的下限
%
%
% 其中，va_upper和va_below 是已知项
% 在单步计算中，已知：
% va_upper
% va_below
% 
% 求出：
% T
% vb_upper
% vb_below
%
% 单独计算完后，逆循环单步中已知：
% T 
% vb
% va_upper
% va_below
%
% 求出：
% va
% vc




% xx_upper 和 xx_below 统称 xx_range
% 单步计算流程：
% STEP 1. va_range -> T_range
% STEP 2. T_range -> T
% STEP 3. va_range, T -> vb_range
% 


% 整体计算流程
% STEP    |  SOURCE
% 1        FOR begin_node : end_node
% 1.1        va_range -> T_range
% 1.2        T_range -> T
% 1.3        va_range, T -> vb_range
% 1        END_FOR
%
% 2        FOR end_node : begin_node
% 2.1        va_range, T, vb -> va, vc
% 2        END_FOR

% STEP 1.1 展开为：
% va_range -> T_range
% 1.1.1    va_upper -> T_below
% 1.1.2    va_below -> T_upper

% STEP 1.2 展开为：
% 1.2.3        BIN_SEARCH T in T_range
%                T = (T_upper + T_below)/2
% 1.2.3          va_range, T -> vb_range
%                FOR next_node : end_node
%                  T = T_upper
% 1.2.3.1          va_range, T -> vb_range
%                  if vb_range contain 0
%                    return true
%                  else if vb_range not exist
%                    return false;
%                  continue;
%                END_FOR                 
%              BIN_END

% STEP 1.3 展开为：
% va_range, T -> vb_range
% 1.3.1 pt, a, j, max_vb -> va_range   (重新计算 va_range)
% 1.3.2 va_upper, T -> vb_upper
% 1.3.3 va_below, T -> vb_below

% STEP 2.1 展开为：
% va_range, T, vb -> va, vc
% ... TBD


% 综上，对于scurve，需要以下函数：
% A s_scurve_cpt_T_upper
%   s_scurve_cpt_T_below
%
% B s_scurve_cpt_vb_upper
%   s_scurve_cpt_vb_below
%
% C s_scurve_cpt_vavc

%% 二、极限情况讨论
%
% 极限情况可以分成四种情况
%
% CASE
% 1    : va -> vc_max -> vb     Tc >= 0
% 2    : va -> vc -> vb         vc > max(va,vb), Tc = 0
% 3    : va -> 0  -> vb         Tc >= 0
% 4    : va -> vc -> vb         vc < min(va,vb), Tc = 0
%
% 以上四种情况，又可以分为四种子情况:
%
% CASE
% A    : va -> vc 无匀加速，vc-> vb 无匀加速
% B    : va -> vc 无匀加速，vc-> vb 有匀加速
% C    : va -> vc 有匀加速，vc-> vb 无匀加速
% D    : va -> vc 有匀加速，vc-> vb 有匀加速
%
%% 计算 vb_upper 与 vb_below
%
% 1. 确定必然可以达到的 vb_max vb_min vc_max vc_min
%
% 2. 确定 vb_min 对应的最大可能的vc：vc_l（vb_max 对应 vc_s）
% 
% 3. 确定 vc_max 对应的最小可能vb：vb_s（vc_min 对应 vb_l）
%
% ------------------------- 求 vb_below -------------------------
%    vc        vb           Ta         Tc       Tb         requires
% 1 vc_max   vb_max       T_va2vcmax T-Ta-Tb T_vcmax2vbmax none
% 2 vc_max   vc_max-a*a/j T_va2vcmax T-Ta-Tb 2*a/j         
% 3 vc_max   vb_s         T_va2vcmax T-Ta-Tb T_vcmax2vbs   none    
% 4 va/-Ta   vc-a*a/j     T-Tb       0       2*a/j         
% 5 va+a*a/j vc\-Tb       2*a/j      0       T-Ta          
% 6 vc_l     vb_min        Ta        Tc         Tb         none
%
% l1 -> l3 -> l6 依次减小
%
% 上式先计算 l1 l3 l6
% l2 若 vc_max-a*a/j > vb_max, 则 l2 = l1
%    若 vc_max-a*a/j < vb_min, 则 l2 = l3
%    否则 按照公式计算
%
% l4 按照公式计算
%
% l5 按照公式计算
% 
% l1 <= pt      时：vb_upper = vb_max
% l2 <= pt < l1 时：有匀速段，Tb段无匀加速，计算Tb
% l3 <= pt < l2 时：有匀速段，Tb段有匀加速，计算Tb
%   T > 4*a*j 时：l4 > l5
%     l4 <= pt < l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l5 <= pt < l4 时：无匀速段，Ta段有匀加速，Tb段有匀加速
%     l6 <= pt < l5 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 4*a*j && T > 2*a*j 时：l5 > l4
%     l5 <= pt < l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l4 <= pt < l5 时：无匀速段，Ta段无匀加速，Tb段无匀加速
%     l6 <= pt < l4 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 2*a*j 时：
%     l6 <= pt < l3 时：无匀速段，Ta段无匀加速，Tb段无匀加速
% pt < l6 时：vb_upper = vb_min
%    
% ------------------------- 求 vb_upper -------------------------
%
%    vc        vb           Ta         Tc       Tb         requires
% 1 vc_min   vb_min       T_va2vcmin T-Ta-Tb T_vcmin2vbmin none
% 2 vc_min   vc_min-a*a/j T_va2vcmin T-Ta-Tb 2*a/j         T, vb
% 3 vc_min   vb_l         T_va2vcmin T-Ta-Tb T_vcmin2vbl   none    
% 4 va/-Ta   vc-a*a/j     T-Tb       0       2*a/j         T > 4*a/j
% 5 va-a*a/j vc\-Tb       2*a/j      0       T-Ta          T > 4*a/j
% 6 vc_s     vb_max       Ta         Tc         Tb         none
%
% l1 -> l3 -> l6 依次增大
%
% 上式先计算 l1 l3 l6
% l2 若 vc_min+a*a/j > vb_max, 则 l2 = l3
%    若 vc_min+a*a/j < vb_min, 则 l2 = l1
%    否则 按照公式计算
%
% l4 按照公式计算
%
% l5 按照公式计算
%
% l1 >= pt      时：vb_upper = vb_min
% l2 >= pt > l1 时：有匀速段，Tb段无匀加速，计算Tb
% l3 >= pt > l2 时：有匀速段，Tb段有匀加速，计算Tb
%   T >  4*a*j 时：l4 < l5
%     l4 >= pt > l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l5 >= pt > l4 时：无匀速段，Ta段有匀加速，Tb段有匀加速
%     l6 >= pt > l5 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 4*a*j && T > 2*a*j 时：l5 < l4
%     l5 >= pt > l3 时：无匀速段，Ta段有匀加速，Tb段无匀加速
%     l4 >= pt > l5 时：无匀速段，Ta段无匀加速，Tb段无匀加速
%     l6 >= pt > l4 时：无匀速段，Ta段无匀加速，Tb段有匀加速
%   T <= 2*a*j 时：
%     l6 >= pt > l3 时：无匀速段，Ta段无匀加速，Tb段无匀加速
% pt > l6 时：vb_upper = vb_min



%
%
% 
%   acc?  vc      vb          Ta    Tc   Tb         vb           requires
% 1  /\ vc_max  vb_max           
% 2  /\ vc_max  vc_max-a*a/j      
% 3  /\ vc_max  vb_s                   
% 4     vc_l    vb_min
%%
% 字母释意：
% 1：加速、有匀速、到达vbmax
% 2：加速、有匀速、到达vbmax - a*a/j
% 3：加速、有匀速、到达0
% 4：加速、有匀速、Tc = 0 vb = 0

% \- 表示减速
% /- 表示加速

%   acc?  vc       Ta        Tc Tb         vb           requires
% 1  /\ vc_max    T_va_vcmax >0 T_vc_vbmax vb_max       T>T_va_vcmax+T_vcmax_vbmax
% 2  /\ vc_max    T_va_vcmax >0 2*a/j      vc_max-a*a/j T>T_va_vcmax+2*a/j,0<vc_max-a*a/j<vb_max  
% 3  /\ vc_max    T_va_vcmax >0 T_0_vcmax  0            T>T_va_vcmax+T_vcmax_0
% 4  /\ vc_max    T_va_vcmax =0 T-Ta       vc_max\-Tb   T>T_va_vcmax
% 5  /\ va/-Ta    T          =0 0          vc           T<T_va_vcmax
% 6  /\ va/-Ta    T-2*a/j    =0 2*a/j      vc-a*a/j     T>4*a/j
% 7  /\ va+a*a/j  2*a/j      =0 T-Ta       >vc-a*a/j    T>4*a/j
% 8  /\ va+a*a/j  2*a/j      =0 T-Ta       <vc-a*a/j    T<4*a/j
% 9  /\ va/-Ta    T-2*a/j    =0 2*a/j      vc-a*a/j     T<4*a/j
% 10  \ va        0          =0 T          va/-T        
% 11  / va        0          =0 T          va\-T        
% 12 \/ va\-Ta    T-2*a/j    =0 2*a/j      vc+a*a/j     T<4*a/j
% 13 \/ va-a*a/j  2*a/j      =0 T-Ta       >vc+a*a/j    T<4*a/j
% 14 \/ va-a*a/j  2*a/j      =0 T-Ta       <vc+a*a/j    T>4*a/j
% 15 \/ va\-Ta    T-2*a/j    =0 2*a/j      vc+a*a/j     T>4*a/j
% 16 \/ va\-Ta    T          =0 0          vc           T<T_va_0
% 17 \/  0        T_va_0     =0 T-Ta       0/-Tb        T>T_va_0
% 18 \/  0        T_va_0     >0 T_0_vbmax  vb_max       T>T_va_0 + T_0_vbmax
% 19 \/  0        T_va_0     >0 2*a/j      a*a/j        T>T_va_0 + 2*a/j
% 20 \/  0        T_va_0     >0 0          0            T>T_va_0
%
% l1  -> l10 依次减小
% l11 -> l20 依次减小


% (vb_max < vc_max-a*a/j时，等同1)

% A : 整个曲线为加速状态（即先加速后减速）
% B : 整个曲线为减速状态（即先减速后加速）
%
% C : 曲线包含匀速段（A时，匀速段为 vc_max, B时为0）
% D : 曲线不包含匀速段
%
% 
%
%
% 分成如下情况：
% 
% 














