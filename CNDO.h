#pragma once
#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include "molecule.h"
#include <cassert>
using namespace std;
using namespace arma;

const double hartree2eV = 27.2114079527;

class CNDO {
    public:
        CNDO(Molecule molecule, mat overlapMatrix);
        Molecule molecule;
        mat overlapMatrix;
        mat gammaMatrix;
        mat hCoreMat;
        mat alphaFockMat;
        mat betaFockMat;
        mat alphaDensityMat;
        mat betaDensityMat;
        mat oldAlphaDensityMat;
        mat oldBetaDensityMat;
        vec totalDensity;
        mat alphaCoeffMat;
        mat betaCoeffMat;
        vec alphaEnergy;
        vec betaEnergy;
        double nuclearRepulsionEnergy;
        double totalEnergy;

        // Map of semi-empirical parameters 
        map<string, map<string, double>> diagCNDOPara;
        map<string, double> offDiagCNDOPara;

        // Map the AO index to the atom it belongs to
        map<int, int> aoIndexToAtom;

        void scfCycle();

    private:
        int calcNumAOs(string chemSym);

        // Functions used to calculate gamma (2e- integral)
        double calc2eIntegral(AO AO1, AO AO2);
        double pg2eIntegral(rowvec center_a, rowvec center_b, double sigmaA, double sigmaB);

        mat calcGammaMatrix();
        vec calcTotalDensity(); 
        mat calcHCoreMat();
        mat calcFockMat(mat densityMat);
        mat calcDensityMat(mat coeffMatA, string type);
        
        double calcNuclearRepulsionEnergy();
        double calcTotalEnergy();
}; 