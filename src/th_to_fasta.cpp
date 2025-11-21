#include <bits/stdc++.h>

using namespace std;

struct Data {
    string instance;
    string column11_sequence;
    double column3_value;
    double column7_value;
};

int main(int argc, const char * argv[]) 
{
    ifstream readsFile, inputFile, lengthFile;
    ofstream outputFile;

    string reads_fname = argv[1];
    string input_fname = argv[2];
    string length_fname = argv[3];
    string output_fname = argv[4];

    readsFile.open(reads_fname.c_str());
    inputFile.open(input_fname.c_str());
    lengthFile.open(length_fname.c_str());
    outputFile.open(output_fname.c_str());

    // Map to store estimated lengths from the 4th file
    map<string, double> estimated_lengths;

    string line, instance;
    
    // Read the estimated lengths from the 4th file
    while (getline(lengthFile, line)) {
        if (line[0] == '>') {
            stringstream ss(line);
            string read_name;
            ss >> read_name; // Extract first token (read name)
            read_name = read_name.substr(1); // Remove '>'
            getline(lengthFile, line); // Read the estimated length (next line)
            double estimated_length = stod(line);
            estimated_lengths[read_name] = estimated_length;
        }
    }

    // Maps to store sequences:
    // - with closest 7th column value to the estimated length
    // - or maximum 3rd column value when no estimated length
    map<string, pair<double, string>> instance_sequences_closest;
    map<string, pair<double, string>> instance_sequences_max3;

    // Read the SRR read headers and store in a vector
    vector<string> headers;
    while (getline(readsFile, line)) {
        if (line[0] == '>') {
            headers.push_back(line);
        }
    }

    // Read the input file and find best candidates
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string temp, column11_sequence;
        double column3_value, column7_value;

        ss >> instance; // Read instance (first column)
        ss >> temp; // Skip 2nd column
        ss >> column3_value; // Read 3rd column

        // Skip to the 7th column
        for (int i = 4; i <= 6; ++i) {
            ss >> temp;
        }
        ss >> column7_value; // Read 7th column

        // Skip to the 11th column
        for (int i = 8; i <= 10; ++i) {
            ss >> temp;
        }
        ss >> column11_sequence; // Read 11th column (sequence)

        if (estimated_lengths.find(instance) != estimated_lengths.end()) {
            double estimated_length = estimated_lengths[instance];
            double current_difference = abs(column7_value - estimated_length);

            if (instance_sequences_closest.find(instance) == instance_sequences_closest.end() || 
                abs(instance_sequences_closest[instance].first - estimated_length) > current_difference) {
                instance_sequences_closest[instance] = {column7_value, column11_sequence};
            }
        } else {
            if (instance_sequences_max3.find(instance) == instance_sequences_max3.end() ||
                instance_sequences_max3[instance].first < column3_value) {
                instance_sequences_max3[instance] = {column3_value, column11_sequence};
            }
        }
    }

    // Write to output file
    for (const string& header : headers) {
        stringstream ss(header);
        string first_token;
        ss >> first_token; // Extract first token (read name)
        string extracted_token = first_token.substr(1); // Remove '>'

        if (instance_sequences_closest.find(extracted_token) != instance_sequences_closest.end()) {
            outputFile << header << endl;
            outputFile << instance_sequences_closest[extracted_token].second << endl;
        } else if (instance_sequences_max3.find(extracted_token) != instance_sequences_max3.end()) {
            outputFile << header << endl;
            outputFile << instance_sequences_max3[extracted_token].second << endl;
        } else {
            outputFile << header << endl;
            outputFile << "norepeat" << endl;
        }
    }

    // Close files
    readsFile.close();
    inputFile.close();
    lengthFile.close();
    outputFile.close();

    return 0;
}
