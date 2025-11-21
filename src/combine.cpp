#include <bits/stdc++.h>
using namespace std;

// Global ordered list of reads
vector<string> read_order;
unordered_set<string> seen_reads;

// Extract read name from a FASTA/cigar header
string extract_readname(const string &line) {
    string s = line.substr(1);     // remove ">"
    size_t pos = s.find(' ');
    if (pos == string::npos) return s;
    return s.substr(0, pos);
}

// Parse cigar value; return -1 if missing or "*"
int extract_cigar_value(const string &line) {
    size_t pos = line.find("cigar_match=");
    if (pos == string::npos) return -1;

    string s = line.substr(pos + 12); // after "cigar_match="

    // Check for "*"
    if (!s.empty() && s[0] == '*') return -1;

    size_t semi = s.find(';');

    // Single value case
    if (semi == string::npos) {
        return stoi(s);
    }

    // Two values: take max
    int a = stoi(s.substr(0, semi));
    int b = stoi(s.substr(semi + 1));
    return max(a, b);
}

// Load cigar file and preserve order of first appearance
unordered_map<string, int> load_cigar(const string &file) {
    unordered_map<string, int> mp;
    mp.reserve(500000);

    ifstream in(file);
    if (!in) {
        cerr << "Error opening cigar file: " << file << endl;
        exit(1);
    }

    string line, current;
    while (getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            current = extract_readname(line);

            // Record order ONCE globally
            if (!seen_reads.count(current)) {
                seen_reads.insert(current);
                read_order.push_back(current);
            }
        }
        else if (line.find("cigar_match=") != string::npos) {
            mp[current] = extract_cigar_value(line);
        }
        else if (line.find("read has *") != string::npos) {
            mp[current] = -1;  // invalid
        }
    }
    return mp;
}

// Load FASTA fully: header + full sequence
unordered_map<string, pair<string,string>> load_fasta(const string &file) {
    unordered_map<string, pair<string,string>> mp;
    mp.reserve(500000);

    ifstream in(file);
    if (!in) {
        cerr << "Error opening FASTA: " << file << endl;
        exit(1);
    }

    string line, current_read, header, seq;

    while (getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current_read.empty()) {
                mp[current_read] = {header, seq};
            }
            header = line;
            current_read = extract_readname(line);
            seq.clear();
        }
        else {
            seq += line;
        }
    }

    if (!current_read.empty()) {
        mp[current_read] = {header, seq};
    }

    return mp;
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        cerr << "Usage: " << argv[0]
             << " A1.cigar A2.cigar A1.fasta A2.fasta output.fasta\n";
        return 1;
    }

    string cigar1_file = argv[1];
    string cigar2_file = argv[2];
    string fasta1_file = argv[3];
    string fasta2_file = argv[4];
    string out_file    = argv[5];

    // Load data (this also fills read_order)
    auto cig1 = load_cigar(cigar1_file);
    auto cig2 = load_cigar(cigar2_file);
    auto fa1  = load_fasta(fasta1_file);
    auto fa2  = load_fasta(fasta2_file);

    ofstream out(out_file);
    if (!out) {
        cerr << "Error opening output file\n";
        return 1;
    }

    // Process in original order
    for (const string &read : read_order) {
        int score1 = cig1.count(read) ? cig1[read] : -1;
        int score2 = cig2.count(read) ? cig2[read] : -1;

        bool useA1 = false;

        if (score1 >= 0 && score2 >= 0) {
            useA1 = (score1 >= score2);
        }
        else if (score1 >= 0) {
            useA1 = true;
        }
        else if (score2 >= 0) {
            useA1 = false;
        }
        else {
            // Both missing or "*"
            continue;
        }

        if (useA1) {
            auto it = fa1.find(read);
            if (it != fa1.end())
                out << it->second.first << "\n" << it->second.second << "\n";
        } else {
            auto it = fa2.find(read);
            if (it != fa2.end())
                out << it->second.first << "\n" << it->second.second << "\n";
        }
    }

    return 0;
}
