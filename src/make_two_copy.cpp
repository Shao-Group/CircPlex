#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_fasta> <output_fasta>" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    ifstream inFile(inputFile);
    ofstream outFile(outputFile);

    if (!inFile || !outFile) {
        cerr << "Error opening input or output file." << endl;
        return 1;
    }

    string line;
    while (getline(inFile, line)) {
        if (line[0] == '>') {
            outFile << line << endl;
            string sequence;
            if (getline(inFile, sequence)) {
                outFile << sequence + sequence << endl;
            }
        }
    }

    inFile.close();
    outFile.close();
    cout << "Duplicated sequences written to " << outputFile << endl;

    return 0;
}
