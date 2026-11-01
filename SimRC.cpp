#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <locale>
using namespace std;

constexpr double E0 = 8.854e-12;
constexpr double e = 1.6e-19;
constexpr double me = 9.109e-31;

const double Cu_n = 8.47e28;
const double Cu_rho0 = 1.68e-8;
const double Cu_alpha = 0.00393;
const double Al_n = 18.1e28;
const double Al_rho0 = 2.65e-8;
const double Al_alpha = 0.00429;
double n;
double rho;
double alpha;
int mat_choice;
double T = 25.0;

//nhập
double L = 1.0;
double Aw = 1e-6;
double V;
double Rm = 1000.0;
double Er = 1.0;
double A_tu = 1e-4;
double d_tu = 1e-3;
double C = 0.0;
double U = 5.0;
int mach_choice; // 1 = mach RC, 2 = test ohm

int main(int argc, char* argv[]) {
    if (argc > 2) {
        cerr << "Cach dung: SimRC.exe [duong_dan_file.csv]\n";
        return 1;
    }
    const filesystem::path outputPath = argc == 2 ? argv[1] : "simrc_data.csv";
    constexpr int N = 400000; 
    constexpr int totalSteps = 8 * N;
    constexpr int consoleInterval = 10000;
    cout << scientific << setprecision(12);

    cout << "Nhap chieu dai day dan (L, m): ";
    cin >> L;
    cout << "Nhap dien tich day dan (Aw, m^2): ";
    cin >> Aw;
    cout << "Nhap dien ap (V, V): ";
    cin >> U;
    cout << "Nhap dien tro (Rm, Ohm): ";
    cin >> Rm;
    if (!cin || !isfinite(L) || !isfinite(Aw) || !isfinite(U)
        || !isfinite(Rm) || L <= 0.0 || Aw <= 0.0 || Rm < 0.0) {
        cerr << "Du lieu khong hop le: can L > 0, Aw > 0, Rm >= 0 va cac so huu han.\n";
        return 1;
    }

    cout << "Chon mach: RC - co tu dien (1) hay thuan tro - khong tu, kiem tra dinh luat Ohm (2): ";
    cin >> mach_choice;
    if (!cin || (mach_choice != 1 && mach_choice != 2)) {
        cerr << "Lua chon mach khong hop le. Chi nhan 1 hoac 2.\n";
        return 1;
    }

    //chọn mode
    int choice = 1; // 1 = nap, 2 = xa
    double Uc0 = 0.0;

    if (mach_choice == 1) {
        cout << "Nhap he so dien moi (Er)[1]: ";
        cin >> Er;
        cout << "Nhap dien tich tu dien (A_tu, m^2): ";
        cin >> A_tu;
        cout << "Nhap khoang cach tu dien (d_tu, m): ";
        cin >> d_tu;
        C = E0 * Er * A_tu / d_tu;
        if (!cin || !isfinite(Er) || !isfinite(A_tu) || !isfinite(d_tu)
            || Er <= 0.0 || A_tu <= 0.0 || d_tu <= 0.0
            || !isfinite(C) || C <= 0.0) {
            cerr << "Du lieu tu dien khong hop le: Er, A_tu, d_tu va C phai duong, huu han.\n";
            return 1;
        }

        cout << "Nap (1) hay xa (2): ";
        if (!(cin >> choice) || (choice != 1 && choice != 2)) {
            cerr << "Lua chon nap/xa khong hop le. Chi nhan 1 hoac 2.\n";
            return 1;
        }
        if (choice == 1) {
            V = U;
            Uc0 = 0.0;
        }
        else {
            V = 0.0;
            Uc0 = U;
        }
    }
    else {
        V = U;
        Uc0 = 0.0;
    }

    //nhập vât liệu dây
    cout << "Chon vat lieu day: Dong ( 1 ), nhom ( 2 )";
    if (!(cin >> mat_choice) || (mat_choice != 1 && mat_choice != 2)) {
        cerr << "Lua chon vat lieu khong hop le. Chi nhan 1 hoac 2.\n";
        return 1;
    }
    if (mat_choice == 1) {
        n = Cu_n;
        rho = Cu_rho0;
        alpha = Cu_alpha;
    }
    else {
        n = Al_n;
        rho = Al_rho0;
        alpha = Al_alpha;
    }

    //nhiệt 
    cout << "Nhap nhiet do (do C, 0-100): ";
    cin >> T;
    double rho_T = rho * (1.0 + alpha * (T - 20.0));
    if (!cin || !isfinite(T) || !isfinite(rho_T) || rho_T <= 0.0) {
        cerr << "Du lieu khong hop le \n";
        return 1;
    }

    //tính biến 
    double tau = (me) / (n * e * e * rho_T);
    double b = -((n * e * e * Aw * Rm) / (me * L)) - (1.0 / tau);
    double Rd = rho_T * L / Aw;
    double dt;
    if (mach_choice == 1) {
        dt = (Rm + Rd) * C / N;
    }
    else {
        dt = (-5.0 / b) / N;
    }
    if (!isfinite(tau) || tau <= 0.0 || !isfinite(b) || b >= 0.0
        || !isfinite(Rd) || Rd <= 0.0 || !isfinite(dt) || dt <= 0.0
        || !isfinite(totalSteps * dt)) {
        cerr << "Thong so tinh toan vuot gioi han so hoc, khong the xuat du lieu.\n";
        return 1;
    }

    double ak = 0.0;
    double vd = 0.0;
    double Qk = Uc0 * C;
    double Uc = Uc0;
    double t_real = 0.0;
    if (!isfinite(Qk)) {
        cerr << "Trang thai ban dau vuot gioi han so hoc, khong the xuat du lieu.\n";
        return 1;
    }

    ofstream dataFile(outputPath);
    if (!dataFile) {
        cerr << "Khong the mo file de ghi: " << outputPath << '\n';
        return 1;
    }
    dataFile.imbue(locale::classic());
    dataFile << scientific << setprecision(numeric_limits<double>::max_digits10);
    dataFile << "step,time_s,Uc_V,vd_m_per_s,I\n";
    // Tai t=0 chua co buoc thoi gian de tinh DeltaQ/dt.
    dataFile << 0 << ',' << t_real << ',' << Uc << ',' << vd << ",\n";

    cout << "t real: " << t_real << " s   ";
    cout << "  Uc: " << Uc << " V   vd: " << vd << " m/s\n";
    cout << "Dang luu toan bo du lieu vao " << outputPath << "...\n";
    // vòng lặp
    for (int step = 1; step <= totalSteps; step++) {
        ak = -(e / (me * L)) * (V - Uc);

        double deltaQ = -n * e * Aw * (-(ak / b) * dt + (vd + ak / b) * expm1(b * dt) / b);
        Qk += deltaQ;
        vd = -(ak / b) + (vd + ak / b) * exp(b * dt);

        if (mach_choice == 1) {
            Uc = Qk / C;
        }
        else {
            Uc = 0.0;
        }

        double I_step = deltaQ / dt;
        t_real = step * dt;
        
        if (!isfinite(deltaQ) || !isfinite(Qk) || !isfinite(Uc)
            || !isfinite(vd) || !isfinite(I_step)
            || !isfinite(t_real)) {
            cerr << "Ket qua khong huu han tai buoc " << step
                << ". File " << outputPath << " chi chua du lieu truoc buoc nay.\n";
            return 1;
        }

        dataFile << step << ',' << t_real << ',' << Uc << ',' << vd << ',' << I_step << '\n';
        if (!dataFile) {
            cerr << "Loi ghi file " << outputPath << " tai buoc " << step
                << ". Du lieu xuat chua day du.\n";
            return 1;
        }

        if (step == 1 || step % consoleInterval == 0 || step == totalSteps) {
            cout << "step: " << step << "  t: " << t_real << " s  Uc: " << Uc
                << " V   vd: " << vd << " m/s  I = DeltaQ/dt: " << I_step << " A\n";
        }
    }

    dataFile.close();
    if (!dataFile) {
        cerr << "Loi hoan tat file " << outputPath << ". Du lieu xuat chua day du.\n";
        return 1;
    }
    error_code pathError;
    const auto fullPath = filesystem::absolute(outputPath, pathError);
    cout << "Da xuat " << totalSteps + 1 << " mau du lieu vao "
        << (pathError ? outputPath : fullPath).string() << '\n';
    return 0;
}
