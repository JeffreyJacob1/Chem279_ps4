#include "CNDO.h"


 // CNDO/2 Hartree Fock method constructor
 
CNDO::CNDO(Molecule molecule, mat overlapMatrix)
    : molecule(molecule), overlapMatrix(overlapMatrix), 
      alphaCoeffMat(arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions)),
      betaCoeffMat(arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions)),
      alphaDensityMat(arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions)),
      betaDensityMat(arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions)),
      totalDensity(arma::zeros(molecule.nBasisFunctions)),
      alphaEnergy(molecule.nBasisFunctions),
      betaEnergy(molecule.nBasisFunctions)
{
    // Initialize semi-empirical parameters
    const vector<string> atoms = {"H", "C", "N", "O", "F"};
    const vector<string> orbitals = {"1s", "2s", "2px", "2py", "2pz"};
    const vector<double> diagValues = {7.176, 14.051, 5.572, 5.572, 5.572, 19.316, 7.275, 7.275, 7.275,
                                       25.390, 9.111, 9.111, 9.111, 32.272, 11.080, 11.080, 11.080};
    const vector<double> offDiagValues = {9., 21., 25., 31., 39.};

    int index = 0;
    for (const string& atom : atoms) {
        if (atom == "H") {
            diagCNDOPara[atom][orbitals[0]] = diagValues[index++];
        } else {
            for (int i = 1; i < orbitals.size(); ++i) {
                diagCNDOPara[atom][orbitals[i]] = diagValues[index++];
            }
        }
        offDiagCNDOPara[atom] = offDiagValues[distance(atoms.begin(), find(atoms.begin(), atoms.end(), atom))];
    }

    // Map the AO index to the atom it belongs to
    index = 0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        int numAOs = calcNumAOs(molecule.atomicSymbols[A]);
        for (int i = 0; i < numAOs; i++) {
            aoIndexToAtom[index++] = A;
        }
    }

    // Calculate matrices
    gammaMatrix = calcGammaMatrix();
    alphaFockMat = calcFockMat(alphaDensityMat);
    betaFockMat = calcFockMat(betaDensityMat);
    hCoreMat = calcHCoreMat();
    nuclearRepulsionEnergy = calcNuclearRepulsionEnergy();
}

/**
 * Calculate the number of AOs for an atom.
 * 
 * input chemSym The atomic symbol of the atom
 * 
 * out int The number of AOs
 */
int CNDO::calcNumAOs(string chemSym) {
    if (chemSym == "H") {
        return 1;
    }
    else {
        return 4;
    }
}

/**
 * Calculate the 2-electron repulsion integrals for two primitive Gaussian shells.
 *  
 * input The first AO, AO2 The second AO
 * 
 * out double The 2-electron repulsion integral
 */
double CNDO::calc2eIntegral(AO atomicOrbital1, AO atomicOrbital2) {
    if (!(accu(atomicOrbital1.lmn) == 0 && accu(atomicOrbital2.lmn) == 0)) {
        cout << "Error: 2e integrals only implemented for s orbitals" << endl;
        return 0.0;
    }
    
    // Compute the product of contraction coefficients and normalization constants
    vec contractionNormProd1 = atomicOrbital1.contraction_coeffs % atomicOrbital1.norm_constants;
    vec contractionNormProd2 = atomicOrbital2.contraction_coeffs % atomicOrbital2.norm_constants;

    int numExponents = atomicOrbital1.exponents.n_elem;
    double twoElectronIntegral = 0.0;

    // Compute the two-electron integral using nested loops
    for (int i = 0; i < numExponents; i++) {
        double exponent1_i = atomicOrbital1.exponents(i);
        double contractionNormProd1_i = contractionNormProd1(i);

        for (int j = 0; j < numExponents; j++) {
            double exponent1_j = atomicOrbital1.exponents(j);
            double contractionNormProd1_j = contractionNormProd1(j);
            double sigma1 = 1.0 / (exponent1_i + exponent1_j);  // eq 3.10

            for (int k = 0; k < numExponents; k++) {
                double exponent2_k = atomicOrbital2.exponents(k);
                double contractionNormProd2_k = contractionNormProd2(k);

                for (int l = 0; l < numExponents; l++) {
                    double exponent2_l = atomicOrbital2.exponents(l);
                    double contractionNormProd2_l = contractionNormProd2(l);
                    double sigma2 = 1.0 / (exponent2_k + exponent2_l);  // eq 3.10

                    double primitiveIntegral = pg2eIntegral(atomicOrbital1.center, atomicOrbital2.center, sigma1, sigma2);  // eq 3.14
                    twoElectronIntegral += contractionNormProd1_i * contractionNormProd1_j * contractionNormProd2_k * contractionNormProd2_l * primitiveIntegral;
                }
            }
        }
    }

    return twoElectronIntegral;
}

/**
 * Calculates the 2-electron integral over primitive Gaussians.
 * 
 * inputs:
 * center_a 
 * center_b 
 * sigmaA 
 * sigmaB 
 * 
 * out: double 2-electron integral
 */
double CNDO::pg2eIntegral(rowvec center_a, rowvec center_b, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V2 = 1.0 / (sigmaA + sigmaB);

    double distance = norm(center_a - center_b, 2);

    if (distance == 0.0) {
        return U * sqrt(2 * V2) * sqrt(2 / M_PI) * hartree2eV;  // eq 3.15
    } 

    double sqrtT = sqrt(V2) * distance;
    double result = U / distance * erf(sqrtT); 
    return result * hartree2eV;
}

/**
 * Calculate the gamma matrix of 2-electron repulsion integrals.
 * 
 * input: molecule 
 * 
 * out: mat The gamma matrix
 */
mat CNDO::calcGammaMatrix() {
    // Make a list of only basis functions with s orbitals
    vector<AO> sBasisFunctionsList;
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        if (molecule.basisFunctionsList[i].AO_type == "1s" || molecule.basisFunctionsList[i].AO_type == "2s") {
            sBasisFunctionsList.push_back(molecule.basisFunctionsList[i]);
        }
    }

    mat gamma_matrix = arma::zeros<mat>(molecule.nAtoms, molecule.nAtoms);

    // Loop over all s orbital basis function combinations
    for (int i = 0; i < sBasisFunctionsList.size(); i++) {
        for (int j = 0; j < sBasisFunctionsList.size(); j++) {
            gamma_matrix(i, j) = calc2eIntegral(sBasisFunctionsList[i], sBasisFunctionsList[j]);
        }
    }
    return gamma_matrix;
}

/**
 * Calculate the nuclear repulsion energy of the molecule.
 * 
 * return double The nuclear repulsion energy
 */
double CNDO::calcNuclearRepulsionEnergy() {
    double nuclearRepulsionEnergy = 0.0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < A; B++) {
            double distance = norm(molecule.coordinates.row(A) - molecule.coordinates.row(B), 2);
            nuclearRepulsionEnergy += molecule.atomValences(A) * molecule.atomValences(B) / distance;
        }
    }
    return nuclearRepulsionEnergy * hartree2eV;
}

   
/**
 * Calculate the total density vector.
 * 
 * return mat The total density vector
 */
vec CNDO::calcTotalDensity() {
    vec totalDensity = arma::zeros(molecule.nAtoms);

    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        // Get atom associated with mu AO
        int A = aoIndexToAtom[mu];
        totalDensity(A) += alphaDensityMat(mu, mu) + betaDensityMat(mu, mu);
    }

    return totalDensity;
}

/**
 * Calculate the CNDO/2 Fock matrix.
 * 
 * param densityMat The density matrix (alpha or beta)
 * 
 * return mat The CNDO/2 Fock matrix
 */
mat CNDO::calcFockMat(mat densityMat) {
    mat fockMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            // Get atoms associated with mu and nu AOs
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double pAA = totalDensity(A);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                fockMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] + \
                                  ((pAA - ZA) - (densityMat(mu, mu) - 0.5)) * gammaAA;
                
                // Update the diagonal elements of the matrix when A != B
                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double pBB = totalDensity(B);
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        fockMat(mu, nu) += (pBB - ZB) * gammaAB;
                    }
                }
            }

            // Calculate the off-diagonal elements of the matrix
            else {
                fockMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                  / 2.0 * overlapMatrix(mu, nu) - (densityMat(mu, nu) * gammaAB);
            }
        }
    }
    return fockMat;
}

/**
 * Calculate the core Hamiltonian matrix.
 * return mat The core Hamiltonian matrix
 */
mat CNDO::calcHCoreMat() {
    mat hCoreMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    
    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                hCoreMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] - (ZA - 0.5) * gammaAA;

                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        hCoreMat(mu, nu) -= ZB * gammaAB;
                    }
                }
            }
            // Calculate the off-diagonal elements of the matrix
            else {
                hCoreMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                   / 2.0 * overlapMatrix(mu, nu);
            }
        }
    }
    return hCoreMat;
}

/**
 * Calculate the density matrix from the coefficient matrix.
 * 
 * param coeffMat The coefficient matrix to use (alpha or beta)
 * param type The type of matrix (alpha or beta)
 * 
 * return mat The density matrix
 */
mat CNDO::calcDensityMat(mat coeffMatA, string type) {
    mat densityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    if (type == "alpha") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.pAlpha; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }

    else if (type == "beta") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.qBeta; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }
    return densityMat;
}

/**
 * Calculate the total energy.
 * 
 * This function calculates the total energy after convergence is reached.
 * 
 * return double The total energy
 */
double CNDO::calcTotalEnergy() {
    double totalEnergy = 0.0;
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            totalEnergy += (alphaDensityMat(mu, nu) * (hCoreMat(mu, nu) + alphaFockMat(mu, nu)) + \
                            betaDensityMat(mu, nu) * (hCoreMat(mu, nu) + betaFockMat(mu, nu)));
        }
    }
    totalEnergy /= 2.0;
    totalEnergy += nuclearRepulsionEnergy;
    return totalEnergy;
}

/**
 * 
 * This function runs the SCF cycle until convergence is reached. The SCF cycle
 * consists of the following steps:
 * 1) Guess the density matrix = 0
 * 2) Calculate the Fock matrix
 * 3) Diagonalize the Fock matrix and obtain the coefficient matrix
 * 4) Calculate the density matrix
 * 5) Check for convergence
 * 6) If not converged, repeat from step 2
 * 7) If converged, calculate the total energy
 */
void CNDO::scfCycle() {
    // 1) Guess the density matrix
    alphaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);

    bool converged = false;
    int scfCycleCount = 0;  

    // 6) If not converged, repeat from step 2
    while (!converged) {
        scfCycleCount++;
        cout << "---------------------------------------------------------" << endl;
        cout << "SCF cycle: " << scfCycleCount << endl << endl;

        // 2) Calculate the Fock matrix
        alphaFockMat = calcFockMat(alphaDensityMat);
        betaFockMat = calcFockMat(betaDensityMat);

        // 3) Diagonalize the Fock matrix and obtain the coefficient matrix
        eig_sym(alphaEnergy, alphaCoeffMat, alphaFockMat);
        eig_sym(betaEnergy, betaCoeffMat, betaFockMat);

        // Save
        oldAlphaDensityMat = alphaDensityMat;
        oldBetaDensityMat = betaDensityMat;

        // 4) Calculate the density matrix
        alphaDensityMat = calcDensityMat(alphaCoeffMat, "alpha");
        betaDensityMat = calcDensityMat(betaCoeffMat, "beta");
        totalDensity = calcTotalDensity();

         // Print matrices
        cout << "Alpha Density Matrix: " << endl;
        cout << alphaDensityMat << endl;
        cout << "Beta Density Matrix: " << endl;
        cout << betaDensityMat << endl;
        cout << "Total Density: " << endl;
        cout << totalDensity << endl;
        cout << "Alpha Fock Matrix: " << endl;
        cout << alphaFockMat << endl;
        cout << "Beta Fock Matrix: " << endl;
        cout << betaFockMat << endl;
        cout << "Alpha Coeff Matrix: " << endl;
        cout << alphaCoeffMat << endl;
        cout << "Beta Coeff Matrix: " << endl;
        cout << betaCoeffMat << endl;

        // 5) Check for convergence (tolerance = 1e-6)
        if (abs(alphaDensityMat - oldAlphaDensityMat).max() < 1e-6 && \
            abs(betaDensityMat - oldBetaDensityMat).max() < 1e-6) {
    
            // 7) If converged, calculate the total energy
            converged = true;
            cout << "*************************************************" << endl;
            cout << "SCF cycle converged after " << scfCycleCount << " iterations!" << endl;
            totalEnergy = calcTotalEnergy();
            cout << "Nuclear Repulsion Energy: " << nuclearRepulsionEnergy << " eV" << endl;
            cout << "Total Energy: " << totalEnergy << " eV" << endl;
        }
    }
}
