#include <iostream>
#include <fstream>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <tuple>

using namespace std;

int extract_unitlength(const string& header) {
    smatch match;
    regex unitlength_re("unitlength=(\\d+)");
    if (regex_search(header, match, unitlength_re)) {
        return stoi(match[1]);
    }
    return -1;
}

int extract_match_from_cigar(const string& cigar) {
    // regex re("([0-9]+)([=M])");
    regex re("([0-9]+)([MI])");
    auto begin = sregex_iterator(cigar.begin(), cigar.end(), re);
    auto end = sregex_iterator();
    int match_total = 0;
    for (auto it = begin; it != end; ++it) {
        int len = stoi((*it)[1]);
        match_total += len;
    }
    return match_total;
}

// ===== NEW FUNCTION =====
// Calculate alignment end from POS and CIGAR string
int calculate_alignment_end(int pos, const string& cigar) {
    regex re("([0-9]+)([MIDNSHP=X])");
    auto begin = sregex_iterator(cigar.begin(), cigar.end(), re);
    auto end = sregex_iterator();
    int ref_len = 0;
    for (auto it = begin; it != end; ++it) {
        int len = stoi((*it)[1]);
        char op = (*it)[2].str()[0];
        if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') {
            ref_len += len;
        }
    }
    return pos + ref_len - 1;
}
// =========================

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " input.fasta input.sam output.log" << endl;
        return 1;
    }

    string fasta_file = argv[1];
    string sam_file = argv[2];
    string log_file = argv[3];

    // Read FASTA: key = read ID, value = {sequence, unitlength, full_header}
    unordered_map<string, tuple<string, int, string>> fasta_map;
    vector<string> read_order;

    ifstream fasta_in(fasta_file);
    string header, seq;

    while (getline(fasta_in, header)) {
        getline(fasta_in, seq);
        if (!header.empty()) {
            int unitlen = extract_unitlength(header);
            string id = header.substr(1);
            id = id.substr(0, id.find(' '));
            fasta_map[id] = make_tuple(seq, unitlen, header);
            read_order.push_back(id);
        }

    }
    fasta_in.close();

    // Read SAM and accumulate match lengths
    unordered_map<string, vector<int>> match_map;
    unordered_map<string, vector<string>> chr_map;
    unordered_map<string, vector<string>> genome_pos_map;
    unordered_map<string, vector<int>> genome_end_map;  // <-- NEW map for end positions

    ifstream sam_in(sam_file);
    string line;
    while (getline(sam_in, line)) {
        if (line.empty() || line[0] == '@') continue;
        vector<string> fields;
        size_t pos = 0, prev = 0;
        while ((pos = line.find('\t', prev)) != string::npos) {
            fields.push_back(line.substr(prev, pos - prev));
            prev = pos + 1;
        }
        fields.push_back(line.substr(prev));

        if (fields.size() < 6) continue;
        string read_id = fields[0];
        string cigar = fields[5];
        string chrm = fields[2];
        string genome_pos = fields[3];

        if (cigar != "*") {
            int match_len = extract_match_from_cigar(cigar);
            match_map[read_id].push_back(match_len);
            chr_map[read_id].push_back(chrm);
            genome_pos_map[read_id].push_back(genome_pos);

            // calculate alignment end and store it
            int pos_start = stoi(genome_pos);
            int pos_end = calculate_alignment_end(pos_start, cigar);
            genome_end_map[read_id].push_back(pos_end);
        }
    }
    sam_in.close();

    // Output in the original FASTA order
    ofstream log_out(log_file);
    for (const string& read_id : read_order) {
        const auto& data = fasta_map[read_id];
        int unitlen = get<1>(data);
        const string& full_header = get<2>(data);

        auto it = match_map.find(read_id);
        if (it == match_map.end()) {
            log_out << full_header << "\nread has * in cigar field\n";
            continue;
        }

        log_out << full_header << "\n";  // print header once

        const vector<int>& matches = match_map[read_id];
        const vector<string>& chrs = chr_map[read_id];
        const vector<string>& genome_positions = genome_pos_map[read_id];
        const vector<int>& genome_ends = genome_end_map[read_id];

        size_t n = matches.size();

        // Join all values with semicolons
        auto join_ints = [](const vector<int>& v) {
            string out;
            for (size_t i = 0; i < v.size(); ++i) {
                if (i > 0) out += ";";
                out += to_string(v[i]);
            }
            return out;
        };
        auto join_strs = [](const vector<string>& v) {
            string out;
            for (size_t i = 0; i < v.size(); ++i) {
                if (i > 0) out += ";";
                out += v[i];
            }
            return out;
        };

        // print compact lines
        log_out << "cigar_match=" << join_ints(matches) << endl;
        log_out << "chr=" << join_strs(chrs) << endl;
        log_out << "pos=" << join_strs(genome_positions) << endl;
        log_out << "end=" << join_ints(genome_ends) << endl;
    }
    log_out.close();

    return 0;
}
