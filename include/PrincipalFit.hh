#pragma once

#include <TH2D.h>
#include <TF1.h>
#include <array>
#include <memory>
// #include <Eigen/Dense>

class PrincipalFit {
public:
    PrincipalFit(int id, TH2D* hist2d);
    ~PrincipalFit() {};
    void SetRegion(double minX_, double maxX_, double minY_, double maxY_) {
        min_x = minX_;
        max_x = maxX_;
        min_y = minY_;
        max_y = maxY_;
        isSet = true;
    };

    void Print();
    void Execution();
    void Draw(int ifunc = 0, const char* option = "same");
    void ShowParameter();
private:
    bool isSet = false;
    std::vector<double> dataX;
    std::vector<double> dataY;
    void GetData();
    int id_;
    double min_x;
    double min_y;
    double max_x;
    double max_y;
    TH2D* h2;
    bool CheckWithinRegion(double x, double y) {
        return min_x <= x && x <= max_x && min_y <= y && y <= max_y;
    };
    void GetCovMatrix();
    void Solve();
    template <class T, class U>
    double GetCovariant(T& input1, U& input2);
    void ObtainEigenVectors();
    std::vector<double> eigenValuesVector;
    std::array<std::vector<double>, 2> eigenVectors;
    std::array<std::shared_ptr<TF1>, 2> funcs;
    void SetParameter();

};

