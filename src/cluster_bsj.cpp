#include <bits/stdc++.h>
using namespace std;

struct BSJ {
    string readname, chr;
    int start, end;
};

// ---- Helper: median of integers ----
int median(vector<int> v) {
    if (v.empty()) return 0;
    sort(v.begin(), v.end());
    int n = v.size();
    return (n % 2) ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
}

// ---- Load a BSJ file ----
vector<BSJ> load_bsj_file(const string &path) {
    ifstream fin(path);
    vector<BSJ> bsjs;
    if (!fin) {
        cerr << "Error: cannot open " << path << endl;
        return bsjs;
    }
    string line;
    getline(fin, line); // skip header
    while (getline(fin, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        BSJ b;
        ss >> b.readname >> b.chr >> b.start >> b.end;
        if (!b.readname.empty()) bsjs.push_back(b);
    }
    return bsjs;
}

// ---- Cluster BSJs by chr/start/end proximity (any match within cluster) ----
vector<vector<BSJ>> cluster_bsj(vector<BSJ> &bsjs, int x) {
    vector<vector<BSJ>> clusters;
    sort(bsjs.begin(), bsjs.end(), [](const BSJ &a, const BSJ &b) {
        if (a.chr != b.chr) return a.chr < b.chr;
        if (a.start != b.start) return a.start < b.start;
        return a.end < b.end;
    });

    for (const auto &b : bsjs) {
        bool added = false;
        for (auto &cl : clusters) {
            // Check if this BSJ matches with any BSJ in the cluster
            for (const auto &existing : cl) {
                if (b.chr == existing.chr &&
                    abs(b.start - existing.start) <= x &&
                    abs(b.end - existing.end) <= x) {
                    cl.push_back(b);
                    added = true;
                    break;
                }
            }
            if (added) break; // no need to check other clusters
        }
        if (!added) {
            // start a new cluster
            clusters.push_back({b});
        }
    }
    return clusters;
}

// ---- Generate output filename ----
string derive_output_name(const string &input) {
    string out = input;
    string suffix = ".bsjs.tsv";
    if (out.size() >= suffix.size() &&
        out.substr(out.size() - suffix.size()) == suffix)
        out = out.substr(0, out.size() - suffix.size());
    out += "_unique.bsjs.tsv";
    return out;
}

// ---- Write clustered output ----
void write_clusters(const vector<vector<BSJ>> &clusters, const string &outfile) {
    ofstream fout(outfile);
    if (!fout) {
        cerr << "Error: cannot write to " << outfile << endl;
        return;
    }
    fout << "readnames\tchr\tstart\tend\tsize\n";
    int written = 0;
    for (const auto &cl : clusters) {
        // if (cl.size() < 2) continue; // Only report clusters with >=2 reads

        vector<int> starts, ends;
        vector<string> names;
        string chr = cl.front().chr;
        for (const auto &b : cl) {
            starts.push_back(b.start);
            ends.push_back(b.end);
            names.push_back(b.readname);
        }

        int cstart = (cl.size() == 1) ? starts[0] : median(starts);
        int cend   = (cl.size() == 1) ? ends[0]   : median(ends);

        fout << names[0];
        for (size_t i = 1; i < names.size(); i++)
            fout << "_" << names[i];
        fout << "\t" << chr << "\t" << cstart << "\t" << cend << "\t" << cl.size() << "\n";
        written++;
    }
    fout.close();
    cerr << "    → " << written << " clusters (size >=2) written to " << outfile << "\n";
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <file1.bsjs.tsv> [file2.bsjs.tsv ...] [x=50]\n";
        return 1;
    }

    int x = 50;
    vector<string> inputs;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg.rfind(".tsv") != string::npos)
            inputs.push_back(arg);
        else
            x = stoi(arg); // last non-.tsv argument = x
    }

    for (auto &path : inputs) {
        cerr << "[*] Processing " << path << " ...\n";
        auto bsjs = load_bsj_file(path);
        auto clusters = cluster_bsj(bsjs, x);
        string outpath = derive_output_name(path);
        write_clusters(clusters, outpath);
    }

    cerr << "[✓] Done.\n";
    return 0;
}
