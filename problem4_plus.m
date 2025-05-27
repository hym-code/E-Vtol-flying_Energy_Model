% 清除工作区和全局变量，确保无旧变量干扰
clear all;
clear global;

% 参数定义
global Ct Cd rou A R m g N c aphar_climb_descent v0_climb num_rotors k1 k2 c4 dt h E_battery v0_cruise aphar_cruise k3 c1 c2 v0_descent;
Ct = 0.03;              % 推力系数
Cd = 0.3;               % 阻力系数
rou = 1.18;             % 空气密度 (kg/m^3)
R = 1;                  % 旋翼半径 (m)
A = 2;                  % 单旋翼面积 (m^2)
m = 1300;               % 飞行器质量 (kg)
g = 9.79;               % 重力加速度 (m/s^2)
N = 2;                  % 单旋翼桨叶数
c = 0.1256;             % 桨叶宽度 (m)
aphar_climb_descent = 0; % 爬升和下降阶段倾角 (度)
aphar_cruise = 5;       % 巡航阶段倾角 (度)
v0_climb = 6;           % 爬升速度 (m/s)
v0_cruise = 160/3.6;    % 巡航速度 (m/s)
v0_descent = -6;        % 下降速度 (m/s，负值表示向下)
h = 200;                % 默认爬升/下降高度 (m)，将被 h_quiet 替换
num_rotors = 8;         % 旋翼数量
E_battery = 100000;     % 电池容量 (Wh)

% 功耗系数
k1 = 0.8824;
k2 = 0.3161;
k3 = 4.7564;
c1 = 0.65;
c2 = 0.487;
c4 = 0.5412;

% 时间步长
dt = 1;                 % 时间步长 (s)

% 噪声模型参数
T_env = 25;             % 环境温度 (°C)
hr = 50;                % 相对湿度 (%)
SPL_day = 55;           % 昼间噪声限制 (dB(A))
SPL_night = 45;         % 夜间噪声限制 (dB(A))
h_min = 200;            % 最小飞行高度 (m)
h_max = 500;            % 最大飞行高度 (m)
h_step = 10;            % 高度步长 (m)

% 对称矩阵（距离，单位：km）
dist_matrix = [
    0.00, 13.16, 26.19, 25.80, 69.68, 51.83, 86.87, 137.96, 150.46;
    13.16, 0.00, 15.99, 17.44, 60.15, 47.67, 84.12, 135.43, 145.78;
    26.19, 15.99, 0.00, 4.19, 44.23, 34.14, 70.72, 121.76, 130.79;
    25.80, 17.44, 4.19, 0.00, 43.90, 30.93, 67.57, 118.75, 128.44;
    69.68, 60.15, 44.23, 43.90, 0.00, 31.33, 48.52, 91.84, 93.91;
    51.83, 47.67, 34.14, 30.93, 31.33, 0.00, 36.64, 87.85, 98.69;
    86.87, 84.12, 70.72, 67.57, 48.52, 36.64, 0.00, 51.33, 65.05;
    137.96, 135.43, 121.76, 118.75, 91.84, 87.85, 51.33, 0.00, 26.85;
    150.46, 145.78, 130.79, 128.44, 93.91, 98.69, 65.05, 26.85, 0.00
];

% 地点名称
locations = {'萧山机场', '钱塘金沙湖', '西湖文化广场', '上城吴山广场', ...
             '临安人民广场', '秦望广场', '桐庐中心广场', '新安江广场', '千岛湖广场'};

% 初始化功耗百分比矩阵
n_locations = length(locations);
power_percent_matrix_day = zeros(n_locations, n_locations);
power_percent_matrix_night = zeros(n_locations, n_locations);

% 计算不扰民高度（昼间和夜间）
[h_quiet_day, n_climb, T_rotor] = calc_quiet_height(SPL_day);
[h_quiet_night, ~, ~] = calc_quiet_height(SPL_night);
fprintf('昼间不扰民高度: %.2f 米\n', h_quiet_day);
fprintf('夜间不扰民高度: %.2f 米\n', h_quiet_night);

% 计算昼间能耗占比
fprintf('\n=== 昼间功耗百分比 ===\n');
t_climb = h_quiet_day / v0_climb;
time_steps_climb = floor(t_climb);
t_descent = h_quiet_day / abs(v0_descent);
time_steps_descent = floor(t_descent);
E_climb = calc_climb_energy(h_quiet_day);
E_descent = calc_descent_energy(h_quiet_day);
for i = 1:n_locations
    for j = 1:n_locations
        if i == j
            power_percent_matrix_day(i, j) = 0;
            continue;
        end
        distance = dist_matrix(i, j);
        t_cruise = (distance * 1000) / v0_cruise;
        time_steps_cruise = floor(t_cruise);
        fprintf('计算 %s 到 %s: distance=%.2f km, time_steps_cruise=%d\n', ...
                locations{i}, locations{j}, distance, time_steps_cruise);
        E_cruise = calc_cruise_energy(time_steps_cruise);
        E_total = E_climb + E_cruise + E_descent;
        power_percent_matrix_day(i, j) = (E_total / E_battery) * 100;
    end
end

% 输出昼间功耗百分比表
fprintf('\n昼间功耗百分比表 (%%):\n');
fprintf('起点/终点\t');
for j = 1:n_locations
    fprintf('%s\t', locations{j});
end
fprintf('\n');
for i = 1:n_locations
    fprintf('%s\t', locations{i});
    for j = 1:n_locations
        fprintf('%.2f\t', power_percent_matrix_day(i, j));
    end
    fprintf('\n');
end

% 计算夜间能耗占比
fprintf('\n=== 夜间功耗百分比 ===\n');
t_climb = h_quiet_night / v0_climb;
time_steps_climb = floor(t_climb);
t_descent = h_quiet_night / abs(v0_descent);
time_steps_descent = floor(t_descent);
E_climb = calc_climb_energy(h_quiet_night);
E_descent = calc_descent_energy(h_quiet_night);
for i = 1:n_locations
    for j = 1:n_locations
        if i == j
            power_percent_matrix_night(i, j) = 0;
            continue;
        end
        distance = dist_matrix(i, j);
        t_cruise = (distance * 1000) / v0_cruise;
        time_steps_cruise = floor(t_cruise);
        fprintf('计算 %s 到 %s: distance=%.2f km, time_steps_cruise=%d\n', ...
                locations{i}, locations{j}, distance, time_steps_cruise);
        E_cruise = calc_cruise_energy(time_steps_cruise);
        E_total = E_climb + E_cruise + E_descent;
        power_percent_matrix_night(i, j) = (E_total / E_battery) * 100;
    end
end

% 输出夜间功耗百分比表
fprintf('\n夜间功耗百分比表 (%%):\n');
fprintf('起点/终点\t');
for j = 1:n_locations
    fprintf('%s\t', locations{j});
end
fprintf('\n');
for i = 1:n_locations
    fprintf('%s\t', locations{i});
    for j = 1:n_locations
        fprintf('%.2f\t', power_percent_matrix_night(i, j));
    end
    fprintf('\n');
end

% 计算不扰民高度函数
function [h_quiet, n_climb, T_rotor] = calc_quiet_height(SPL_target)
    global Ct Cd rou A R m g N c aphar_climb_descent v0_climb num_rotors;
    % 噪声模型参数
    T_env = 25;         % 环境温度 (°C)
    hr = 50;            % 相对湿度 (%)
    h_min = 200;        % 最小飞行高度 (m)
    h_max = 500;        % 最大飞行高度 (m)
    h_step = 10;        % 高度步长 (m)
    h_values = h_min:h_step:h_max;
    
    % 计算爬升阶段转速和推力
    n_guess = 2000;
    fun = @(n) (num_rotors * Ct * rou * A * (2 * pi * n/60 * R)^2) - ...
               (m*g + calc_resistance(n, Ct, Cd, rou, A, R, v0_climb, N, c, aphar_climb_descent));
    options = optimoptions('fsolve', 'Display', 'none', 'Algorithm', 'levenberg-marquardt', ...
                          'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                          'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
    try
        n_climb = fsolve(fun, n_guess, options);
        n_climb = max(1000, min(4000, n_climb));
        if ~isscalar(n_climb)
            error('calc_quiet_height: fsolve 返回非标量值，n = %s', mat2str(n_climb));
        end
    catch e
        fprintf('calc_quiet_height: fsolve 失败，错误: %s\n', e.message);
        h_quiet = h_min;
        n_climb = 2000;
        T_rotor = 0;
        return;
    end
    
    w = 2 * pi * n_climb / 60;
    T_temp = Ct * rou * A * (w * R)^2;
    T_rotor = (m * g + calc_resistance(n_climb, Ct, Cd, rou, A, R, v0_climb, N, c, aphar_climb_descent)) / num_rotors;
    
    % 搜索满足噪声约束的最小高度
    h_quiet = h_max; % 默认最大高度
    for h = h_values
        SPL = fw_h_noise_with_absorption(h, v0_climb, T_env, hr, n_climb, T_rotor);
        if SPL <= SPL_target
            h_quiet = h;
            fprintf('高度 %.2f 米, SPL = %.2f dB, 满足目标 %s dB\n', h, SPL, num2str(SPL_target));
            break;
        end
    end
    if h_quiet == h_max
        fprintf('警告: 在 %d-%d 米范围内未找到满足 SPL <= %s dB 的高度，采用默认高度 %.2f 米\n', ...
                h_min, h_max, num2str(SPL_target), h_quiet);
    end
end

% 噪声模型函数（修改自第一个代码）
function [SPL] = fw_h_noise_with_absorption(h, u_inf, T, hr, n, T_rotor)
    % 物理常数
    T_K = T + 273.15; % 转换为绝对温度（K）
    c0 = 331.3 * sqrt(T_K / 273.15); % 音速 (m/s)
    rho0 = 1.18; % 空气密度 (kg/m^3)
    p_ref = 2e-5; % 参考声压 (Pa)

    % 旋翼参数
    R = 1; % 旋翼半径 (m)
    Omega = n * 2 * pi / 60; % 动态转速 (rad/s)
    B = 2; % 刀片数
    P_surface = T_rotor / (pi * R^2); % 表面压力基于推力 (Pa)
    mu = 1; % 谐波序数
    n_mode = 0; % 非稳态模式

    % 计算马赫数
    M_t = Omega * R / c0;
    J = u_inf / (Omega * R);

    % 进距比影响函数
    G_J = -175 - 8*J + 185*J^2;

    % 大气吸收系数（调整以避免负SPL）
    f = 1000; % 假设频率
    alpha = 0.005 * (T_K / 293.15)^(-1/2) * (50 / hr); % 减小吸收系数

    % 观察点（距离h，径向）
    r = h;
    theta = pi/2;

    % 权重系数（从第一个代码调整）
    k1 = 1.45;
    k2 = 1.9;
    k3 = 2.1;
    k4 = 1.2;

    % FW-H模型：载荷噪声
    if h >= 200 && h < 350
        p = k4 * ((1/(4*pi*r)) * B * P_surface * (M_t^(2 + abs(mu*B + n_mode))) * G_J);
    else
        p = (1/(4*pi*r)) * B * P_surface * (M_t^(2 + abs(mu*B + n_mode))) * G_J;
    end

    % 考虑大气吸收
    SPL = k1 * 20 * log10(p / p_ref) - k2 * alpha * r - k3 * 20 * log10(2 * pi * h);
end

% 爬升阶段能量计算函数（接受动态高度）
function E_climb = calc_climb_energy(h_quiet)
    global Ct Cd rou A R m g N c aphar_climb_descent v0_climb num_rotors k1 k2 c4 dt;
    % 检查全局变量并打印非标量或未定义变量
    vars = {'Ct', 'Cd', 'rou', 'A', 'R', 'm', 'g', 'N', 'c', 'aphar_climb_descent', ...
            'v0_climb', 'num_rotors', 'k1', 'k2', 'c4', 'dt'};
    for i = 1:length(vars)
        if ~exist(vars{i}, 'var') || isempty(eval(vars{i}))
            error('calc_climb_energy: 全局变量 %s 未定义或为空', vars{i});
        end
        if ~isscalar(eval(vars{i}))
            error('calc_climb_energy: 全局变量 %s 必须是标量，当前值: %s', ...
                  vars{i}, mat2str(eval(vars{i})));
        end
    end
    time_steps_climb = floor(h_quiet / v0_climb);
    if time_steps_climb <= 0
        E_climb = 0;
        return;
    end
    P_total_climb = zeros(1, time_steps_climb);
    
    for i = 1:time_steps_climb
        n_guess = 2000; % 初始猜测值
        fun = @(n) (num_rotors * Ct * rou * A * (2 * pi * n/60 * R)^2) - ...
                   (m*g + calc_resistance(n, Ct, Cd, rou, A, R, v0_climb, N, c, aphar_climb_descent));
        options = optimoptions('fsolve', 'Display', 'none', 'Algorithm', 'levenberg-marquardt', ...
                              'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                              'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
        try
            n = fsolve(fun, n_guess, options);
            n = max(1000, min(4000, n));
            if ~isscalar(n)
                error('calc_climb_energy: fsolve 返回非标量值，n = %s', mat2str(n));
            end
        catch e
            fprintf('calc_climb_energy: fsolve 失败，错误: %s\n', e.message);
            E_climb = 0;
            return;
        end
        
        w = 2 * pi * n / 60;
        T_temp = Ct * rou * A * (w * R)^2;
        v1 = sqrt(T_temp / (2 * rou * A));
        v2 = v0_climb + v1;
        f1 = 0.5 * rou * A * (v2^2 - v1^2);
        f2 = 0.5 * Cd * rou * A * (w * R)^2;
        f3 = 0.5 * rou * v0_climb^2 * A * Cd;
        f = f1 + f2 + f3;
        T = (m * g + f) / num_rotors;
        
        P1 = k1 * T * (v0_climb / 2 + ((v0_climb / 2)^2 + T / (k2)^2)^(1/2));
        P3 = c4 * (v0_climb)^3;
        P = (P1 + P3) / (60 * 0.8);
        P_total_climb(i) = num_rotors * P;
    end
    
    E_climb = sum(P_total_climb * dt) / 3600;
end

% 巡航阶段能量计算函数
function E_cruise = calc_cruise_energy(time_steps_cruise)
    global Ct Cd rou A R m g N c aphar_cruise v0_cruise num_rotors k1 k2 k3 c1 c2 c4 dt;
    % 检查全局变量并打印非标量或未定义变量
    vars = {'Ct', 'Cd', 'rou', 'A', 'R', 'm', 'g', 'N', 'c', 'aphar_cruise', ...
            'v0_cruise', 'num_rotors', 'k1', 'k2', 'k3', 'c1', 'c2', 'c4', 'dt'};
    for i = 1:length(vars)
        if ~exist(vars{i}, 'var') || isempty(eval(vars{i}))
            error('calc_cruise_energy: 全局变量 %s 未定义或为空', vars{i});
        end
        if ~isscalar(eval(vars{i}))
            error('calc_cruise_energy: 全局变量 %s 必须是标量，当前值: %s', ...
                  vars{i}, mat2str(eval(vars{i})));
        end
    end
    if abs(cosd(aphar_cruise)) < 1e-6
        error('calc_cruise_energy: cosd(aphar_cruise) 接近零，可能导致数值错误');
    end
    if time_steps_cruise <= 0
        E_cruise = 0;
        return;
    end
    P_total_cruise = zeros(1, time_steps_cruise);
    T = zeros(1, time_steps_cruise);
    for i = 1:time_steps_cruise
        n_guess = 2000;
        fun = @(n) (num_rotors * Ct * rou * A * (2 * pi * n/60 * R)^2 * cosd(aphar_cruise)) - ...
                   (m*g + calc_resistance(n, Ct, Cd, rou, A, R, v0_cruise, N, c, aphar_cruise));
        options = optimoptions('fsolve', 'Display', 'none', 'Algorithm', 'levenberg-marquardt', ...
                              'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                              'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
        try
            n = fsolve(fun, n_guess, options);
            n = max(1000, min(4000, n));
            if ~isscalar(n)
                error('calc_cruise_energy: fsolve 返回非标量值，n = %s', mat2str(n));
            end
        catch e
            fprintf('calc_cruise_energy: fsolve 失败，错误: %s\n', e.message);
            E_cruise = 0;
            return;
        end
        
        w = 2 * pi * n / 60;
        resistance = calc_resistance(n, Ct, Cd, rou, A, R, v0_cruise, N, c, aphar_cruise);
        T(i) = (m * g + resistance) / (num_rotors * cosd(aphar_cruise));
        
        v1 = sqrt(T(i) / (2 * rou * A));
        v2 = v0_cruise + v1;
        P1 = k1 * T(i) * (v0_cruise / 2 + ((v0_cruise / 2)^2 + T(i) / (k2)^2)^(1/2));
        mu = (v0_cruise * cosd(aphar_cruise)) / (w * R);
        P2 = (c * Cd * rou * R^4 * w^3) * (1 + mu^2);
        P3 = c4 * (v0_cruise)^3;
        P = ((c1 + c2) * T(i)^(3/2) + P1 + P3 + k3 * P2) / (3600 * 0.8);
        P_total_cruise(i) = num_rotors * P;
    end
    
    E_cruise = sum(P_total_cruise * dt) / 3600;
end

% 下降阶段能量计算函数（接受动态高度）
function E_descent = calc_descent_energy(h_quiet)
    global Ct Cd rou A R m g N c aphar_climb_descent v0_descent num_rotors k1 k2 c4 dt;
    % 检查全局变量并打印非标量或未定义变量
    vars = {'Ct', 'Cd', 'rou', 'A', 'R', 'm', 'g', 'N', 'c', 'aphar_climb_descent', ...
            'v0_descent', 'num_rotors', 'k1', 'k2', 'c4', 'dt'};
    for i = 1:length(vars)
        if ~exist(vars{i}, 'var') || isempty(eval(vars{i}))
            error('calc_descent_energy: 全局变量 %s 未定义或为空', vars{i});
        end
        if ~isscalar(eval(vars{i}))
            error('calc_descent_energy: 全局变量 %s 必须是标量，当前值: %s', ...
                  vars{i}, mat2str(eval(vars{i})));
        end
    end
    time_steps_descent = floor(h_quiet / abs(v0_descent));
    if time_steps_descent <= 0
        E_descent = 0;
        return;
    end
    P_total_descent = zeros(1, time_steps_descent);
    local_num_rotors = num_rotors + 4; % 使用局部变量，避免修改全局 num_rotors
    
    for i = 1:time_steps_descent
        n_guess = 2000;
        fun = @(n) (local_num_rotors * Ct * rou * A * (2 * pi * n/60 * R)^2) - ...
                   (m*g + calc_resistance(n, Ct, Cd, rou, A, R, v0_descent, N, c, aphar_climb_descent));
        options = optimoptions('fsolve', 'Display', 'none', 'Algorithm', 'levenberg-marquardt', ...
                              'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8, ...
                              'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);
        try
            n = fsolve(fun, n_guess, options);
            n = max(1000, min(4000, n));
            if ~isscalar(n)
                error('calc_descent_energy: fsolve 返回非标量值，n = %s', mat2str(n));
            end
        catch e
            fprintf('calc_descent_energy: fsolve 失败，错误: %s\n', e.message);
            E_descent = 0;
            return;
        end
        
        w = 2 * pi * n / 60;
        T_temp = Ct * rou * A * (w * R)^2;
        v1 = sqrt(T_temp / (2 * rou * A));
        v2 = v0_descent + v1;
        f1 = 0.5 * rou * A * (v2^2 - v1^2);
        f2 = 0.5 * Cd * rou * A * (w * R)^2;
        f3 = 0.5 * rou * v0_descent^2 * A * Cd;
        f = f1 + f2 + f3;
        T = (m * g + f) / local_num_rotors;
        
        P1 = k1 * T * (abs(v0_descent) / 2 + ((abs(v0_descent) / 2)^2 + T / (k2)^2)^(1/2));
        P3 = c4 * (abs(v0_descent))^3;
        P = (P1 + P3) / (60 * 2);
        P_total_descent(i) = local_num_rotors * P;
    end
    
    E_descent = sum(P_total_descent * dt) / 3600;
end

% 阻力计算函数
function f = calc_resistance(n, Ct, Cd, rou, A, R, v0, N, c, aphar)
    if ~isscalar(n) || ~isscalar(v0) || ~isscalar(aphar)
        fprintf('calc_resistance: 非标量输入检测到:\n');
        fprintf('n = %s\n', mat2str(n));
        fprintf('v0 = %s\n', mat2str(v0));
        fprintf('aphar = %s\n', mat2str(aphar));
        error('calc_resistance: 输入参数 n, v0, aphar 必须是标量');
    end
    w = 2 * pi * n / 60;
    T = Ct * rou * A * (w * R)^2;
    v1 = sqrt(T / (2 * rou * A));
    v2 = v0 + v1;
    f1 = 0.5 * rou * A * (v2^2 - v1^2);
    f2 = 0.5 * Cd * rou * A * (w * R)^2;
    f3 = 0.5 * rou * v0^2 * A * Cd;
    f = f1 + f2 + f3;
    if ~isfinite(f)
        fprintf('calc_resistance: 输出 f 不是有限值，f = %s\n', mat2str(f));
        error('calc_resistance: 输出 f 不是有限值');
    end
end
