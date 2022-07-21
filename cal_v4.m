clc
clear all
%% Parameters
phi = 0.06;  % 孔隙率，1
L   = 6.5*10^(-3); % 样品长度，m
R   = 1.1*10^(-3); % 样品截面半径，m
v_p = pi * R^2 * L * phi; % 样品孔隙体积，m^3
v_u = 690.72 * 10^(-9); % 上游体积，m^3
v_d = 698.14 * 10^(-9); % 下游体积，m^3
a   = v_p / v_u / 8; % 样品孔隙体积/上游体积 / 8 ，1
b   = v_p / v_d / 8; % 样品孔隙体积/下游体积 / 8 ，1 
mu  = 1.79 * 10^(-5); % 氮气动力粘性系数，26.4℃，1bar
nu  = 1.59 * 10^(-5); % 氮气运动粘性系数，26.4℃，1bar
%% 样品编号 39#, phi=0.06
t_iner_tot_temp  = 0; % 记录 inertial solution 总用时
t_late_tot_temp  = 0; % 记录 late-time solution 总用时 
for ii = 1 : 1
    exp1 = xlsread( strcat(num2str(ii),'.xlsx')); % 读入文件
    data_delta_t = 10;  % 数据记录时间间隔，s
    % Inertial solution, no slip  
    t_seq_iner   = [];
    dpD_seq_iner = []; 
    begin_iner   = 1;
    for jj = begin_iner : length(exp1)
        temp_iner = log( exp((exp1(jj,3) / ( exp1(begin_iner,3) - exp1(begin_iner,1) ))/16) - exp((exp1(jj,1) / ( exp1(begin_iner,3) - exp1(begin_iner,1) ))/16));
        temp = ( exp1(jj,3) - exp1(jj,1) ) / ( exp1(1,3) - exp1(1,1) );
        t_seq_iner(jj - begin_iner + 1,1)   = data_delta_t + (jj - begin_iner) * data_delta_t; 
        dpD_seq_iner(jj - begin_iner + 1,1) = temp_iner; 
        if temp < 0.95
            index_iner_end = jj;
            break
        end
    end
    syms t_iner
    f = fittype(' - A * t_iner + B ', 'independent', 't_iner', 'coefficients', {'A','B'} );
    [cfun_iner_sin, output_iner_sin] = fit(t_seq_iner, dpD_seq_iner, f);
    
    t_iner_tot(ii,1) = t_iner_tot_temp + max(t_seq_iner);
    rsquare_iner(ii,1) = output_iner_sin.rsquare;  % 拟合线性相关系数
    alpha(ii,1) = cfun_iner_sin.A;
    sita1 = fzero(@(sita1) (sita1 .^2 - a * b ) .* tan(sita1)  - (a + b) * sita1 ,[0.0000001 0.1]);
    P_u0_iner(ii,1) = exp1(begin_iner,3) * 10^5; 
    P_d0_iner(ii,1) = exp1(begin_iner,1) * 10^5; 
    kappa_iner(ii,1) =  1 / 8 * phi * mu  * L^2 * alpha(ii,1) / sita1^2 / P_u0_iner(ii,1);
    
    pr1(ii,1) = kappa_iner(ii,1) * 8 / phi * P_u0_iner(ii,1) / mu / nu;
    pr2(ii,1) = (exp1(1,3) - exp1(1,1)) / exp1(1,3);
    pr3(ii,1) = pr1(ii,1) * pr2(ii,1);
    
    % Late time solution, no slip
%     for jj = 1 : length(exp1)
%         temp_late(jj,1)          = (exp1(jj,3) - exp1(jj,1)) / exp1(jj,3);
%     end
%     [row1, col1] = find ( abs ( temp_late - 0.15 ) < 0.001);
%     [row2, col2] = find ( abs ( temp_late - 0.10 ) < 0.001);
%     
%     dpD_seq_late = [];
%     t_seq_late   = [];
%     P_u0_late(ii,1)    = exp1(min(row1),3) * 10^5;
%     for jj = min(row1) : max(row2)
%         dpD_seq_late( (jj - min(row1) + 1),1) = log ( (exp1(jj,3) - exp1(jj,1)) / ( exp1(min(row1),3) - exp1(min(row1),1)) );
%     end
%     t_seq_late(:,1)    = [0 : data_delta_t : data_delta_t * (length(dpD_seq_late)-1) ];
%     syms t_late
%     f = fittype(' C * t_late + D ', 'independent', 't_late', 'coefficients', {'C','D'} );
%     [cfun, output] = fit(t_seq_late, dpD_seq_late, f);
%     t_late_tot(ii,1) = t_late_tot_temp  + max(t_seq_late);
%     mathbb_B(ii,1) = cfun.C;
%     sita1 = fzero(@(sita1) (sita1 .^2 - a * b ) .* tan(sita1)  - (a + b) * sita1 ,[0.0000001 0.1]);
%     kappa_late(ii,1) = - mathbb_B(ii,1) * mu * phi * L^2 / 8 / sita1^2 / P_u0_late(ii);
%     
%     sita1_trad = fzero(@(sita1_trad) (sita1_trad .^2 - 64 * a * b ) .* tan(sita1_trad)  - 8 * (a + b) * sita1_trad ,[0.0000001 0.1]);
%     kappa_trad(ii,1) = - mathbb_B(ii,1) * mu * phi * L^2 / sita1_trad^2 / P_u0_late(ii);
end

% syms P_rev
% f = fittype(' E * P_rev + F ', 'independent', 'P_rev', 'coefficients', {'E','F'} );
% [cfun_late, output_late] = fit( 1./P_u0_late(3:6,1), kappa_late(3:6,1)*1E19, f); 
% kappa_late_infty(1:length(kappa_late),1) = cfun_late.F;

kappa_iner_of1(1:ii,1) = spline(pr3, kappa_iner, 1.0);
%% output
% T = table(P_u0_iner, P_d0_iner, P_u0_late, kappa_trad, kappa_late, kappa_late_infty, kappa_iner, kappa_iner_of1, pr1, pr2, pr3, rsquare_iner, t_iner_tot, t_late_tot);
% writetable(T,'output.xlsx');
% type output.xlsx;
T = table(P_u0_iner, kappa_iner, kappa_iner_of1, pr1, pr2, pr3, rsquare_iner, t_iner_tot);
writetable(T,'output.xlsx');
type output.xlsx;