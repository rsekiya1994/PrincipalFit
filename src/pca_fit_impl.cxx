#include "../macro/pca_fit.C"
#include <TROOT.h>
#include <memory>
#include <TH2D.h>
#include <TSystem.h>
#include "PrincipalFit.hh"

std::shared_ptr<PrincipalFit> pca_ptr;

void pca_fit(const char* histname, double min_x, double max_x, double min_y, double max_y) {
    pca_ptr.reset();
    TH2D* h2 = (TH2D*)gROOT->FindObject(histname);
    pca_ptr = std::make_shared<PrincipalFit>(1, h2);
    pca_ptr->SetRegion(min_x, max_x, min_y, max_y);
    pca_ptr->Execution();
    pca_ptr->Draw(0);
    pca_ptr->Draw(1);
    pca_ptr->ShowParameter();
}