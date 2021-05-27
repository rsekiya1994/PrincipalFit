#include "PrincipalFit.hh"
#include <iostream>
#include <numeric>
#include "Eigen/Dense"
#include <TString.h>

struct SolverStr {
    Eigen::Matrix2d covMatrix;
    Eigen::EigenSolver<Eigen::Matrix<double, 2, 2> > solver;
} solv;

PrincipalFit::PrincipalFit(int id, TH2D* hist2d) :id_(id) {
    h2 = hist2d;
}

void PrincipalFit::Execution() {
    if (!isSet) {
        std::cout << "Error : Range is not set " << std::endl;
        return;
    }
    GetData();
    GetCovMatrix();
    Solve();
    ObtainEigenVectors();
    funcs[0] = std::make_shared<TF1>(Form("func%d", id_), "pol1", min_x, max_x);
    funcs[1] = std::make_shared<TF1>(Form("func%d", -id_), "pol1", min_x, max_x);
    SetParameter();
}

void PrincipalFit::GetData() {
    int nBinX = h2->GetXaxis()->GetNbins();
    double widthX = h2->GetXaxis()->GetBinWidth(0);
    double minX = h2->GetXaxis()->GetXmin();
    double maxX = h2->GetXaxis()->GetXmax();

    int nBinY = h2->GetYaxis()->GetNbins();
    double widthY = h2->GetYaxis()->GetBinWidth(0);
    double minY = h2->GetYaxis()->GetXmin();
    double maxY = h2->GetYaxis()->GetXmax();
    int nEntries = h2->GetEntries();
    // std::cout << "==== INFO ====\n";
    // std::cout << "-- num bin X = " << nBinX << "\n";
    // std::cout << "-- min X     = " << minX << "\n";
    // std::cout << "-- max X     = " << maxX << "\n";
    // std::cout << "-- width X   = " << widthX << "\n";
    // std::cout << "-- num bin Y = " << nBinY << "\n";
    // std::cout << "-- min Y     = " << minY << "\n";
    // std::cout << "-- max Y     = " << maxY << "\n";
    // std::cout << "-- width Y   = " << widthY << "\n";
    // std::cout << "nEntries " << nEntries << "\n";
    // std::cout << std::endl;
    dataX.reserve(nEntries);
    dataY.reserve(nEntries);

    for (int i = 0; i < nBinX; ++i) {
        for (int j = 0; j < nBinY; ++j) {
            double nowX = minX + widthX * i;
            double nowY = minY + widthY * j;
            auto isWithin = CheckWithinRegion(nowX, nowY);
            if (!isWithin) continue;
            int contents = static_cast<int>(h2->GetBinContent(i, j));
            // std::cout << nowX << " " << nowY << std::endl;
            for (int k = 0; k < contents; ++k) {
                dataX.emplace_back(nowX);
                dataY.emplace_back(nowY);
            }
        }// for <j>
    }// for <i>
}

void PrincipalFit::GetCovMatrix() {
    solv.covMatrix(0, 0) = GetCovariant(dataX, dataX);
    solv.covMatrix(1, 0) = solv.covMatrix(0, 1) = GetCovariant(dataX, dataY);
    solv.covMatrix(1, 1) = GetCovariant(dataY, dataY);
}

template <class T, class U>
double PrincipalFit::GetCovariant(T& input1, U& input2) {
    auto size = input1.size();
    auto cov = std::inner_product(input1.begin(), input1.end(), input2.begin(), 0.0) / (double)(size);
    auto ave1 = std::accumulate(input1.begin(), input1.end(), 0.0,
                                [](decltype(*input1.begin()) acc, decltype(*input1.begin()) i) {return acc + i;});
    auto ave2 = std::accumulate(input2.begin(), input2.end(), 0.0,
                                [](decltype(*input2.begin()) acc, decltype(*input2.begin()) i) {return acc + i;});

    ave1 /= (double)(size);
    ave2 /= (double)(size);
    return cov - ave1 * ave2;
}

void PrincipalFit::Print() {
    std::cout << "-- Covariant Matrix \n";
    std::cout << solv.covMatrix << std::endl;
    std::cout << "-- Eigen Vectors \n";
    std::cout << solv.solver.eigenvectors() << "\n";
    std::cout << "-- Eigen Values \n";
    std::cout << solv.solver.eigenvalues() << "\n";
    std::cout << std::endl;
}

void PrincipalFit::Solve() {
    solv.solver.compute(solv.covMatrix);
}

void PrincipalFit::ObtainEigenVectors() {
    Eigen::MatrixXcd eigenVec(2, 2);
    // Eigen::MatrixXcd eigenValues(2, 2);
    // std::cout << "-- 1" << std::endl;
    eigenVec = solv.solver.eigenvectors();
    auto eigenValues = solv.solver.eigenvalues();
    // auto eigenValues1 = solv.solver.eigenvalues();
    // std::cout << "-- 2" << std::endl;
    // eigenValuesVector = { eigenValues(0, 0).real(), eigenValues(1,1).real() };
    // std::cout << eigenValues << std::endl;
    eigenValuesVector = { eigenValues(0).real(), eigenValues(1).real() };
    // std::cout << "-- 3" << std::endl;
    eigenVectors[0] = { eigenVec(0, 0).real(), eigenVec(1, 0).real() };
    // std::cout << "-- 4" << std::endl;
    eigenVectors[1] = { eigenVec(0, 1).real(), eigenVec(1, 1).real() };
}


void PrincipalFit::SetParameter() {
    // y - meanY = p0 + p1 * (x - meanX)
    // y = p0 - meanX * p1 + meanY + p1 * x;
    double p0[2];
    double p1[2];
    double meanX = std::accumulate(dataX.begin(), dataX.end(), 0.0) / dataX.size();
    double meanY = std::accumulate(dataY.begin(), dataY.end(), 0.0) / dataY.size();
    if (eigenValuesVector[0] > eigenValuesVector[1]) {
        p1[0] = eigenVectors[0][1] / eigenVectors[0][0];
        p1[1] = eigenVectors[1][1] / eigenVectors[1][0];
    } else {
        p1[0] = eigenVectors[1][1] / eigenVectors[1][0];
        p1[1] = eigenVectors[0][1] / eigenVectors[0][0];
    }
    p0[0] = meanY - p1[0] * meanX;
    p0[1] = meanY - p1[1] * meanX;
    funcs[0]->SetParameter(0, p0[0]);
    funcs[0]->SetParameter(1, p1[0]);
    funcs[1]->SetParameter(0, p0[1]);
    funcs[1]->SetParameter(1, p1[1]);
}


void PrincipalFit::Draw(int ifunc, const char* option) {
    funcs[ifunc]->Draw(option);
}


void PrincipalFit::ShowParameter() {
    std::sort(eigenValuesVector.rbegin(), eigenValuesVector.rend());
    std::cout << "============ Fit Parameter Info =========== \n";
    std::cout << "------ Func 1 : \n";
    std::cout << "-- p0  = " << funcs[0]->GetParameter(0) << "\n";
    std::cout << "-- p1  = " << funcs[0]->GetParameter(1) << "\n";
    std::cout << "-- Std Dev. = " << std::sqrt(eigenValuesVector[0]) << "\n";
    std::cout << "------ Func 2 : \n";
    std::cout << "-- p0  = " << funcs[1]->GetParameter(0) << "\n";
    std::cout << "-- p1  = " << funcs[1]->GetParameter(1) << std::endl;
    std::cout << "-- Std Dev. = " << std::sqrt(eigenValuesVector[1]) << "\n";
}