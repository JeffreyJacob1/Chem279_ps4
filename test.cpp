#include <iostream>
#include "molecule.h"
#include "CNDO.h"
#include <chrono>
#include <cassert>
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Incorrect file input, ./test Inputs/file.txt << endl";
        return 1;
    }

    // Time the program
    auto start_time = chrono::high_resolution_clock::now();

    const char* filename = argv[1];
    Molecule molecule(filename);
    cout << "Molecule Info: " << endl;
    molecule.printMoleculeInfo();

    cout << "----------------------------------------" << endl;

    // Test all matrix functions before SCF cycle
    cout << "Initial Matrix Calculations: " << endl << endl;
    mat S = calcOverlapMatrix(molecule);
    cout << "Overlap Matrix: " << endl;
    cout << S << endl;
    
    CNDO cndo(molecule, S);
    cout << "Gamma Matrix: " << endl;
    cout << cndo.gammaMatrix << endl;

    cout << "H Core Matrix: " << endl;
    cout << cndo.hCoreMat << endl;

    cout << "Alpha Fock Matrix: " << endl;
    cout << cndo.alphaFockMat << endl;

    cout << "Beta Fock Matrix: " << endl;
    cout << cndo.betaFockMat << endl;

    // Test SCF cycle
    cndo.scfCycle();

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "Time taken: " << duration.count() << " microseconds" << endl;

    return 0;
}