#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <functional>
#include <thread>
#include <mutex>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace std;
using namespace Eigen;

// 全局常量参数
const double Ct = 0.03;              // 推力系数
const double Cd = 0.3;               // 阻力系数
const double rou = 1.18;             // 空气密度 (kg/m^3)
const double R = 1.0;                // 旋翼半径 (m)
const double A = 2.0;                // 单旋翼面积 (m^2)
const double m = 1300.0;             // 飞行器质量 (kg)
const double g = 9.79;               // 重力加速度 (m/s^2)
const int N = 2;                     // 单旋翼桨叶数
const double c = 0.1256;             // 桨叶宽度 (m)
const double aphar_climb_descent = 0.0; // 爬升和下降阶段倾角 (度)
const double aphar_cruise = 5.0;     // 巡航阶段倾角 (度)
const double v0_climb = 6.0;         // 爬升速度 (m/s)
const double v0_cruise = 160.0 / 3.6; // 巡航速度 (m/s)
const double v0_descent = -6.0;      // 下降速度 (m/s，负值表示向下)
const double h = 200.0;              // 默认爬升/下降高度 (m)
const int num_rotors = 8;            // 旋翼数量
const double E_battery = 100000.0;   // 电池容量 (Wh)

// 功耗系数
const double k1 = 0.8824;
const double k2 = 0.3161;
const double k3 = 4.7564;
const double c1 = 0.65;
const double c2 = 0.487;
const double c4 = 0.5412;

// 时间步长
const double dt = 1.0;               // 时间步长 (s)

// 噪声模型参数
const double T_env = 25.0;           // 环境温度 (°C)
const double hr = 50.0;              // 相对湿度 (%)
const double SPL_day = 55.0;         // 昼间噪声限制 (dB(A))
const double SPL_night = 45.0;       // 夜间噪声限制 (dB(A))
const double h_min = 200.0;          // 最小飞行高度 (m)
const double h_max = 500.0;          // 最大飞行高度 (m)
const double h_step = 10.0;          // 高度步长 (m)

// 距离矩阵 (km)
const vector<vector<double>> dist_matrix = {
    {0.00, 13.16, 26.19, 25.80, 69.68, 51.83, 86.87, 137.96, 150.46},
    {13.16, 0.00, 15.99, 17.44, 60.15, 47.67, 84.12, 135.43, 145.78},
    {26.19, 15.99, 0.00, 4.19, 44.23, 34.14, 70.72, 121.76, 130.79},
    {25.80, 17.44, 4.19, 0.00, 43.90, 30.93, 67.57, 118.75, 128.44},
    {69.68, 60.15, 44.23, 43.90, 0.00, 31.33, 48.52, 91.84, 93.91},
    {51.83, 47.67, 34.14, 30.93, 31.33, 0.00, 36.64, 87.85, 98.69},
    {86.87, 84.12, 70.72, 67.57, 48.52, 36.64, 0.00, 51.33, 65.05},
    {137.96, 135.43, 121.76, 118.75, 91.84, 87.85, 51.33, 0.00, 26.85},
    {150.46, 145.78, 130.79, 128.44, 93.91, 98.69, 65.05, 26.85, 0.00}
};

// 地点名称
const vector<string> locations = {
    "萧山机场", "钱塘金沙湖", "西湖文化广场", "上城吴山广场",
    "临安人民广场", "秦望广场", "桐庐中心广场", "新安江广场", "千岛湖广场"
};

// 互斥锁用于线程安全输出
mutex cout_mutex;

// 计算阻力函数（内联）
inline double calc_resistance(const double n, const double Ct, const double Cd, const double rou, const double A,
                             const double R, const double v0, const int N, const double c, const double aphar) {
    if (!isfinite(n) || !isfinite(v0) || !isfinite(aphar)) {
        throw runtime_error("calc_resistance: 检测到非标量或非有限输入");
    }
    double w = 2.0 * M_PI * n / 60.0;
    double T = Ct * rou * A * pow(w * R, 2);
    double v1 = sqrt(T / (2.0 * rou * A));
    double v2 = v0 + v1;
    double f1 = 0.5 * rou * A * (pow(v2, 2) - pow(v1, 2));
    double f2 = 0.5 * Cd * rou * A * pow(w * R, 2);
    double f3 = 0.5 * rou * pow(v0, 2) * A * Cd;
    double f = f1 + f2 + f3;
    if (!isfinite(f)) {
        throw runtime_error("calc_resistance: 输出 f 非有限值");
    }
    return f;
}

// 计算噪声函数（内联）
inline double fw_h_noise_with_absorption(const double h, const double u_inf, const double T, const double hr,
                                        const double n, const double T_rotor) {
    double T_K = T + 273.15; // 转换为绝对温度（K）
    double c0 = 331.3 * sqrt(T_K / 273.15); // 音速
    double rho0 = 1.18; // 空气密度
    double p_ref = 2e-5; // 参考声压
    double R = 1.0; // 旋翼半径
    double Omega = n * 2.0 * M_PI / 60.0; // 角速度
    int B = 2; // 桨叶数
    double P_surface = T_rotor / (M_PI * pow(R, 2)); // 表面压力
    int mu = 1; // 谐波序数
    int n_mode = 0; // 非稳态模式

    // 计算马赫数
    double M_t = Omega * R / c0;
    double J = u_inf / (Omega * R);

    // 进距比影响函数
    double G_J = -175.0 - 8.0 * J + 185.0 * pow(J, 2);

    // 大气吸收系数
    double f = 1000.0; // 假设频率
    double alpha = 0.005 * pow(T_K / 293.15, -0.5) * (50.0 / hr);

    // 观察点
    double r = h;
    double theta = M_PI / 2.0;

    // 权重系数
    double k1 = 1.45;
    double k2 = 1.9;
    double k3 = 2.1;
    double k4 = 1.2;

    // FW-H模型：载荷噪声
    double p;
    if (h >= 200.0 && h < 350.0) {
        p = k4 * ((1.0 / (4.0 * M_PI * r)) * B * P_surface * pow(M_t, 2 + abs(mu * B + n_mode)) * G_J);
    } else {
        p = (1.0 / (4.0 * M_PI * r)) * B * P_surface * pow(M_t, 2 + abs(mu * B + n_mode)) * G_J;
    }

    // 计算带大气吸收的SPL
    return k1 * 20.0 * log10(p / p_ref) - k2 * alpha * r - k3 * 20.0 * log10(2.0 * M_PI * h);
}

// 非线性求解器仿真MATLAB的fsolve
struct Functor {
    const double Ct, rou, A, R, m, g, v0;
    const int N, num_rotors;
    const double c, aphar;
    Functor(const double Ct, const double rou, const double A, const double R, const double m, const double g,
            const double v0, const int N, const double c, const double aphar, const int num_rotors)
        : Ct(Ct), rou(rou), A(A), R(R), m(m), g(g), v0(v0), N(N), c(c), aphar(aphar), num_rotors(num_rotors) {}
    int operator()(const VectorXd &x, VectorXd &fvec) const {
        double n = x(0);
        fvec(0) = (num_rotors * Ct * rou * A * pow(2.0 * M_PI * n / 60.0 * R, 2)) -
                  (m * g + calc_resistance(n, Ct, Cd, rou, A, R, v0, N, c, aphar));
        return 0;
    }
    int df(const VectorXd &x, MatrixXd &fjac) const {
        double n = x(0);
        double h = 1e-8;
        VectorXd xh = x;
        xh(0) += h;
        VectorXd fvec1(1), fvec2(1);
        (*this)(x, fvec1);
        (*this)(xh, fvec2);
        fjac(0, 0) = (fvec2(0) - fvec1(0)) / h;
        return 0;
    }
    int inputs() const { return 1; }
    int values() const { return 1; }
};

// 计算不扰民高度
inline void calc_quiet_height(const double SPL_target, double &h_quiet, double &n_climb, double &T_rotor) {
    vector<double> h_values;
    for (double h = h_min; h <= h_max; h += h_step) {
        h_values.push_back(h);
    }

    // 求解转速
    double n_guess = 2000.0;
    Functor functor(Ct, rou, A, R, m, g, v0_climb, N, c, aphar_climb_descent, num_rotors);
    LevenbergMarquardt<Functor> lm(functor);
    VectorXd x(1);
    x(0) = n_guess;
    lm.minimize(x);
    n_climb = max(1000.0, min(4000.0, x(0)));

    double w = 2.0 * M_PI * n_climb / 60.0;
    double T_temp = Ct * rou * A * pow(w * R, 2);
    T_rotor = (m * g + calc_resistance(n_climb, Ct, Cd, rou, A, R, v0_climb, N, c, aphar_climb_descent)) / num_rotors;

    // 搜索满足噪声约束的最小高度
    h_quiet = h_max;
    for (double h : h_values) {
        double SPL = fw_h_noise_with_absorption(h, v0_climb, T_env, hr, n_climb, T_rotor);
        if (SPL <= SPL_target) {
            h_quiet = h;
            {
                lock_guard<mutex> lock(cout_mutex);
                cout << "高度 " << fixed << setprecision(2) << h << " 米, SPL = " << SPL
                     << " dB, 满足目标 " << SPL_target << " dB" << endl;
            }
            break;
        }
    }
    if (h_quiet == h_max) {
        lock_guard<mutex> lock(cout_mutex);
        cout << "警告: 在 " << h_min << "-" << h_max
             << " 米范围内未找到满足 SPL <= " << SPL_target << " dB 的高度，使用默认高度 " << h_quiet << " 米" << endl;
    }
}

// 计算爬升阶段能耗
inline double calc_climb_energy(const double h_quiet) {
    int time_steps_climb = floor(h_quiet / v0_climb);
    if (time_steps_climb <= 0) return 0.0;

    vector<double> P_total_climb(time_steps_climb, 0.0);
    for (int i = 0; i < time_steps_climb; ++i) {
        double n_guess = 2000.0;
        Functor functor(Ct, rou, A, R, m, g, v0_climb, N, c, aphar_climb_descent, num_rotors);
        LevenbergMarquardt<Functor> lm(functor);
        VectorXd x(1);
        x(0) = n_guess;
        lm.minimize(x);
        double n = max(1000.0, min(4000.0, x(0)));

        double w = 2.0 * M_PI * n / 60.0;
        double T_temp = Ct * rou * A * pow(w * R, 2);
        double v1 = sqrt(T_temp / (2.0 * rou * A));
        double v2 = v0_climb + v1;
        double f1 = 0.5 * rou * A * (pow(v2, 2) - pow(v1, 2));
        double f2 = 0.5 * Cd * rou * A * pow(w * R, 2);
        double f3 = 0.5 * rou * pow(v0_climb, 2) * A * Cd;
        double f = f1 + f2 + f3;
        double T = (m * g + f) / num_rotors;

        double P1 = k1 * T * (v0_climb / 2.0 + sqrt(pow(v0_climb / 2.0, 2) + T / pow(k2, 2)));
        double P3 = c4 * pow(v0_climb, 3);
        double P = (P1 + P3) / (60.0 * 0.8);
        P_total_climb[i] = num_rotors * P;
    }

    double E_climb = 0.0;
    for (double P : P_total_climb) {
        E_climb += P * dt;
    }
    return E_climb / 3600.0;
}

// 计算巡航阶段能耗
inline double calc_cruise_energy(const int time_steps_cruise) {
    if (abs(cos(aphar_cruise * M_PI / 180.0)) < 1e-6) {
        throw runtime_error("calc_cruise_energy: cos(aphar_cruise) 接近零");
    }
    if (time_steps_cruise <= 0) return 0.0;

    vector<double> P_total_cruise(time_steps_cruise, 0.0);
    vector<double> T(time_steps_cruise, 0.0);
    for (int i = 0; i < time_steps_cruise; ++i) {
        double n_guess = 2000.0;
        Functor functor(Ct, rou, A, R, m, g, v0_cruise, N, c, aphar_cruise, num_rotors);
        LevenbergMarquardt<Functor> lm(functor);
        VectorXd x(1);
        x(0) = n_guess;
        lm.minimize(x);
        double n = max(1000.0, min(4000.0, x(0)));

        double w = 2.0 * M_PI * n / 60.0;
        double resistance = calc_resistance(n, Ct, Cd, rou, A, R, v0_cruise, N, c, aphar_cruise);
        T[i] = (m * g + resistance) / (num_rotors * cos(aphar_cruise * M_PI / 180.0));

        double v1 = sqrt(T[i] / (2.0 * rou * A));
        double v2 = v0_cruise + v1;
        double P1 = k1 * T[i] * (v0_cruise / 2.0 + sqrt(pow(v0_cruise / 2.0, 2) + T[i] / pow(k2, 2)));
        double mu = (v0_cruise * cos(aphar_cruise * M_PI / 180.0)) / (w * R);
        double P2 = (c * Cd * rou * pow(R, 4) * pow(w, 3)) * (1.0 + pow(mu, 2));
        double P3 = c4 * pow(v0_cruise, 3);
        double P = ((c1 + c2) * pow(T[i], 1.5) + P1 + P3 + k3 * P2) / (3600.0 * 0.8);
        P_total_cruise[i] = num_rotors * P;
    }

    double E_cruise = 0.0;
    for (double P : P_total_cruise) {
        E_cruise += P * dt;
    }
    return E_cruise / 3600.0;
}

// 计算下降阶段能耗
inline double calc_descent_energy(const double h_quiet) {
    int time_steps_descent = floor(h_quiet / abs(v0_descent));
    if (time_steps_descent <= 0) return 0.0;

    int local_num_rotors = num_rotors + 4;
    vector<double> P_total_descent(time_steps_descent, 0.0);
    for (int i = 0; i < time_steps_descent; ++i) {
        double n_guess = 2000.0;
        Functor functor(Ct, rou, A, R, m, g, v0_descent, N, c, aphar_climb_descent, local_num_rotors);
        LevenbergMarquardt<Functor> lm(functor);
        VectorXd x(1);
        x(0) = n_guess;
        lm.minimize(x);
        double n = max(1000.0, min(4000.0, x(0)));

        double w = 2.0 * M_PI * n / 60.0;
        double T_temp = Ct * rou * A * pow(w * R, 2);
        double v1 = sqrt(T_temp / (2.0 * rou * A));
        double v2 = v0_descent + v1;
        double f1 = 0.5 * rou * A * (pow(v2, 2) - pow(v1, 2));
        double f2 = 0.5 * Cd * rou * A * pow(w * R, 2);
        double f3 = 0.5 * rou * pow(v0_descent, 2) * A * Cd;
        double f = f1 + f2 + f3;
        double T = (m * g + f) / local_num_rotors;

        double P1 = k1 * T * (abs(v0_descent) / 2.0 + sqrt(pow(abs(v0_descent) / 2.0, 2) + T / pow(k2, 2)));
        double P3 = c4 * pow(abs(v0_descent), 3);
        double P = (P1 + P3) / (60.0 * 2.0);
        P_total_descent[i] = local_num_rotors * P;
    }

    double E_descent = 0.0;
    for (double P : P_total_descent) {
        E_descent += P * dt;
    }
    return E_descent / 3600.0;
}

// 多线程计算功耗矩阵的函数
void compute_power_matrix(const int i, vector<vector<double>> &power_matrix, const double h_quiet,
                         const double E_climb, const double E_descent, const bool is_day) {
    for (size_t j = 0; j < locations.size(); ++j) {
        if (i == j) {
            power_matrix[i][j] = 0.0;
            continue;
        }
        double distance = dist_matrix[i][j];
        int time_steps_cruise = floor((distance * 1000.0) / v0_cruise);
        {
            lock_guard<mutex> lock(cout_mutex);
            cout << "计算 " << locations[i] << " 到 " << locations[j]
                 << ": 距离=" << fixed << setprecision(2) << distance
                 << " 公里, 巡航时间步长=" << time_steps_cruise << endl;
        }
        double E_cruise = calc_cruise_energy(time_steps_cruise);
        double E_total = E_climb + E_cruise + E_descent;
        power_matrix[i][j] = (E_total / E_battery) * 100.0;
    }
}

int main() {
    int n_locations = locations.size();
    vector<vector<double>> power_percent_matrix_day(n_locations, vector<double>(n_locations, 0.0));
    vector<vector<double>> power_percent_matrix_night(n_locations, vector<double>(n_locations, 0.0));

    // 计算不扰民高度
    double h_quiet_day, n_climb_day, T_rotor_day;
    double h_quiet_night, n_climb_night, T_rotor_night;
    calc_quiet_height(SPL_day, h_quiet_day, n_climb_day, T_rotor_day);
    calc_quiet_height(SPL_night, h_quiet_night, n_climb_night, T_rotor_night);
    cout << "昼间不扰民高度: " << fixed << setprecision(2) << h_quiet_day << " 米" << endl;
    cout << "夜间不扰民高度: " << h_quiet_night << " 米" << endl;

    // 昼间功耗计算
    cout << "\n=== 昼间功耗百分比 ===\n";
    int time_steps_climb_day = floor(h_quiet_day / v0_climb);
    int time_steps_descent_day = floor(h_quiet_day / abs(v0_descent));
    double E_climb_day = calc_climb_energy(h_quiet_day);
    double E_descent_day = calc_descent_energy(h_quiet_day);

    vector<thread> threads_day;
    for (int i = 0; i < n_locations; ++i) {
        threads_day.emplace_back(compute_power_matrix, i, ref(power_percent_matrix_day),
                                 h_quiet_day, E_climb_day, E_descent_day, true);
    }
    for (auto &t : threads_day) {
        t.join();
    }

    // 输出昼间功耗百分比表
    cout << "\n昼间功耗百分比表 (%):\n";
    cout << "起点/终点\t";
    for (const auto &loc : locations) {
        cout << loc << "\t";
    }
    cout << endl;
    for (int i = 0; i < n_locations; ++i) {
        cout << locations[i] << "\t";
        for (int j = 0; j < n_locations; ++j) {
            cout << fixed << setprecision(2) << power_percent_matrix_day[i][j] << "\t";
        }
        cout << endl;
    }

    // 夜间功耗计算
    cout << "\n=== 夜间功耗百分比 ===\n";
    int time_steps_climb_night = floor(h_quiet_night / v0_climb);
    int time_steps_descent_night = floor(h_quiet_night / abs(v0_descent));
    double E_climb_night = calc_climb_energy(h_quiet_night);
    double E_descent_night = calc_descent_energy(h_quiet_night);

    vector<thread> threads_night;
    for (int i = 0; i < n_locations; ++i) {
        threads_night.emplace_back(compute_power_matrix, i, ref(power_percent_matrix_night),
                                  h_quiet_night, E_climb_night, E_descent_night, false);
    }
    for (auto &t : threads_night) {
        t.join();
    }

    // 输出夜间功耗百分比表
    cout << "\n夜间功耗百分比表 (%):\n";
    cout << "起点/终点\t";
    for (const auto &loc : locations) {
        cout << loc << "\t";
    }
    cout << endl;
    for (int i = 0; i < n_locations; ++i) {
        cout << locations[i] << "\t";
        for (int j = 0; j < n_locations; ++j) {
            cout << fixed << setprecision(2) << power_percent_matrix_night[i][j] << "\t";
        }
        cout << endl;
    }

    return 0;
}