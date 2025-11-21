#include <bits/stdc++.h>
using namespace std;

// Reverse complement of a DNA k-mer
string revcomp(const string &s) {
    string rc;
    rc.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': rc += 'T'; break;
            case 'T': rc += 'A'; break;
            case 'C': rc += 'G'; break;
            case 'G': rc += 'C'; break;
            default:  rc += 'N'; break;
        }
    }
    return rc;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <kmer_count_file> <log_file>" << endl;
        return 1;
    }

    string infile = argv[1];
    string logfile = argv[2];

    ifstream fin(infile);
    ofstream log(logfile, ios::app); // append mode

    if (!fin.is_open()) {
        cerr << "Error: cannot open " << infile << endl;
        return 1;
    }
    if (!log.is_open()) {
        cerr << "Error: cannot open log file " << logfile << endl;
        return 1;
    }

    unordered_map<string, int> kmer_count;
    kmer_count.reserve(1000000);

    string kmer;
    int count;
    long long total_kmer_count = 0;

    // ---- Read kmers ----
    while (fin >> kmer >> count) {
        if (count >= 2) {
            kmer_count[kmer] = count;
            total_kmer_count += count;
        }
    }
    fin.close();

    if (total_kmer_count == 0) {
        log << "[INFO] No kmers with count >=2. Marked as simple.\n";
        cout << "simple" << endl;
        return 0;
    }

    long long rev_kmer_count = 0;
    unordered_set<string> visited;
    visited.reserve(kmer_count.size());

    log << "===== Reverse Complement Pair Analysis =====\n";
    log << "File: " << infile << "\n";
    log << "Total kmers (count>=2): " << kmer_count.size() << "\n\n";

    for (const auto &p : kmer_count) {
        const string &k1 = p.first;
        if (visited.count(k1)) continue;

        string rc = revcomp(k1);
        auto it = kmer_count.find(rc);
        if (it != kmer_count.end()) {
            int c1 = p.second, c2 = it->second;
            int minc = min(c1, c2);
            rev_kmer_count += minc;
            log << k1 << " (" << c1 << ")\t<->\t"
                << rc << " (" << c2 << ")"
                << "  => min=" << minc << "\n";
            visited.insert(rc);
        }
        visited.insert(k1);
    }

    double ratio = (double)rev_kmer_count / (double)total_kmer_count;
    const double THRESHOLD = 0.05;

    log << "\nrev_kmer_count = " << rev_kmer_count
        << ", total_kmer_count = " << total_kmer_count
        << ", ratio = " << fixed << setprecision(3) << ratio
        << ", threshold = " << THRESHOLD << "\n";

    string decision = (ratio >= THRESHOLD) ? "complex" : "simple";
    log << "Decision: " << decision << "\n";
    log << "============================================\n\n";

    cout << decision << endl;
    return 0;
}
