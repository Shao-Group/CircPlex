#include <bits/stdc++.h>
using namespace std;

// Global ordered list of reads
vector<string> read_order;
unordered_set<string> seen_reads;

// Extract read name from a FASTA/cigar header
string extract_readname(const string &line) {
    string s = line.substr(1); // remove ">"
    size_t pos = s.find(' ');
    if (pos == string::npos) return s;
    return s.substr(0, pos);
}

// Split a semicolon-separated string into ints
vector<int> split_ints(const string &s) {
    vector<int> v;
    stringstream ss(s);
    string token;
    while (getline(ss, token, ';')) {
        if (token == "*") v.push_back(-1);
        else v.push_back(stoi(token));
    }
    return v;
}

// Split a semicolon-separated string (like chr) into strings
vector<string> split_strs(const string &s) {
    vector<string> v;
    stringstream ss(s);
    string token;
    while (getline(ss, token, ';')) v.push_back(token);
    return v;
}

// Cigar info for one read
struct CigarInfo {
    vector<int> cigar_match;
    vector<string> chr;
    vector<int> start;
    vector<int> end;
};

// Parse cigar file
unordered_map<string, CigarInfo> load_cigar(const string &file) {
    unordered_map<string, CigarInfo> mp;
    ifstream in(file);
    if (!in) {
        cerr << "Error opening cigar file: " << file << endl;
        exit(1);
    }

    string line, current;
    CigarInfo info;

    while (getline(in, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!current.empty()) mp[current] = info; // save previous
            current = extract_readname(line);
            if (!seen_reads.count(current)) {
                seen_reads.insert(current);
                read_order.push_back(current);
            }
            info = CigarInfo(); // reset
        } else if (line.find("cigar_match=") != string::npos) {
            string s = line.substr(line.find('=')+1);
            info.cigar_match = split_ints(s);
        } else if (line.find("chr=") != string::npos) {
            string s = line.substr(line.find('=')+1);
            info.chr = split_strs(s);
        } else if (line.find("pos=") != string::npos) {
            string s = line.substr(line.find('=')+1);
            info.start = split_ints(s);
        } else if (line.find("end=") != string::npos) {
            string s = line.substr(line.find('=')+1);
            info.end = split_ints(s);
        } else if (line.find("read has *") != string::npos) {
            info.cigar_match = {-1}; // mark as invalid
        }
    }
    if (!current.empty()) mp[current] = info;
    return mp;
}

// Load FASTA fully: header + full sequence
unordered_map<string, pair<string,string>> load_fasta(const string &file) {
    unordered_map<string, pair<string,string>> mp;
    ifstream in(file);
    if (!in) {
        cerr << "Error opening FASTA: " << file << endl;
        exit(1);
    }

    string line, current_read, header, seq;
    while (getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_read.empty()) mp[current_read] = {header, seq};
            header = line;
            current_read = extract_readname(line);
            seq.clear();
        } else seq += line;
    }
    if (!current_read.empty()) mp[current_read] = {header, seq};
    return mp;
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        cerr << "Usage: " << argv[0] 
             << " A1.cigar A2.cigar A1.fasta A2.fasta output.fasta output_bsj.tsv\n";
        return 1;
    }

    string cigar1_file = argv[1];
    string cigar2_file = argv[2];
    string fasta1_file = argv[3];
    string fasta2_file = argv[4];
    string out_file    = argv[5];
    string out_bsj     = argv[6];

    auto cig1 = load_cigar(cigar1_file);
    auto cig2 = load_cigar(cigar2_file);
    auto fa1  = load_fasta(fasta1_file);
    auto fa2  = load_fasta(fasta2_file);

    ofstream out(out_file);
    if (!out) { cerr << "Error opening output fasta\n"; return 1; }

    ofstream bsj(out_bsj);
    if (!bsj) { cerr << "Error opening output BSJ file\n"; return 1; }
    bsj << "readname\tsource\tchr\tstart\tend\n";

    // Process reads in order
    for (const string &read : read_order) {
        CigarInfo c1 = cig1.count(read) ? cig1[read] : CigarInfo{{-1},{},{},{}};
        CigarInfo c2 = cig2.count(read) ? cig2[read] : CigarInfo{{-1},{},{},{}};

        int max1 = c1.cigar_match.empty() ? -1 : *max_element(c1.cigar_match.begin(), c1.cigar_match.end());
        int max2 = c2.cigar_match.empty() ? -1 : *max_element(c2.cigar_match.begin(), c2.cigar_match.end());

        bool useA1 = false;
        int idx = -1;

        if (max1 >= 0 && max2 >= 0) {
            if (max1 >= max2) { useA1 = true; idx = distance(c1.cigar_match.begin(), max_element(c1.cigar_match.begin(), c1.cigar_match.end())); }
            else             { useA1 = false; idx = distance(c2.cigar_match.begin(), max_element(c2.cigar_match.begin(), c2.cigar_match.end())); }
        } else if (max1 >= 0) { useA1 = true; idx = distance(c1.cigar_match.begin(), max_element(c1.cigar_match.begin(), c1.cigar_match.end())); }
        else if (max2 >= 0) { useA1 = false; idx = distance(c2.cigar_match.begin(), max_element(c2.cigar_match.begin(), c2.cigar_match.end())); }
        else continue; // skip invalid

        // Write fasta output
        if (useA1) {
            auto it = fa1.find(read);
            if (it != fa1.end())
                out << it->second.first << "\n" << it->second.second << "\n";
        } else {
            auto it = fa2.find(read);
            if (it != fa2.end())
                out << it->second.first << "\n" << it->second.second << "\n";
        }

        // Write BSJ output
        const CigarInfo &chosen = useA1 ? c1 : c2;
        string chr = (idx >=0 && idx < chosen.chr.size()) ? chosen.chr[idx] : "NA";
        int start  = (idx >=0 && idx < chosen.start.size()) ? chosen.start[idx] : -1;
        int end    = (idx >=0 && idx < chosen.end.size()) ? chosen.end[idx] : -1;

        bsj << read << "\tCircPlex\t" << chr << "\t" << start << "\t" << end << "\n";
    }

    return 0;
}
