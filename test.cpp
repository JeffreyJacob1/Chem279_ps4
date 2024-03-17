#include <iostream>
#include "molecule.h"
#include "CNDO.h"
#include <chrono>
#include <cassert>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Incorrect file input. Usage: ./test Inputs/file.txt" << endl;
        return 1;
    }

    // Start timing the program
    auto startTime = chrono::high_resolution_clock::now();

    // Read the input file
    const char* inputFile = argv[1];
    Molecule molecule(inputFile);

    // Print molecule information
    cout << "Molecule Information:" << endl;
    molecule.printMoleculeInfo();
    cout << "----------------------------------------" << endl;

    // Perform initial matrix calculations
    cout << "Initial Matrix Calculations:" << endl << endl;

    // Calculate the overlap matrix
    mat overlapMatrix = calcOverlapMatrix(molecule);
    cout << "Overlap Matrix:" << endl;
    cout << overlapMatrix << endl;

    // Create CNDO object
    CNDO cndoObject(molecule, overlapMatrix);

    // Print matrices
    cout << "Gamma Matrix:" << endl;
    cout << cndoObject.gammaMatrix << endl;

    cout << "H Core Matrix:" << endl;
    cout << cndoObject.hCoreMat << endl;

    cout << "Alpha Fock Matrix:" << endl;
    cout << cndoObject.alphaFockMat << endl;

    cout << "Beta Fock Matrix:" << endl;
    cout << cndoObject.betaFockMat << endl;

    // Perform SCF cycle
    cndoObject.scfCycle();

    // End timing the program
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(endTime - startTime);

    // Print the execution time
    cout << "Execution Time: " << duration.count() << " microseconds" << endl;

    return 0;
}