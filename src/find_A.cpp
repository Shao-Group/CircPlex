#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <climits>
#include <iomanip>

using namespace std;

#define vvi vector <vi>
#define ALIGNMENT_THRESHOLD 40

string reverse_complement(const string& sequence) {
    std::string complemented_sequence = "";
    for (int i = sequence.size() - 1; i >= 0; i--) {
        char base = sequence[i];
        char complement_base;
        if (base == 'A')
            complement_base = 'T';
        else if (base == 'T')
            complement_base = 'A';
        else if (base == 'C')
            complement_base = 'G';
        else if (base == 'G')
            complement_base = 'C';
        else
            complement_base = base;
        complemented_sequence += complement_base;
    }
    return complemented_sequence;
}

pair<int, int> findClosestPair(const vector<pair<int, int>>& pairs) {
    pair<int, int> closestPair;
    int minDiff = numeric_limits<int>::max();
    for (const auto& pair : pairs) {
        int diff = abs(pair.first - pair.second);
        if (diff < minDiff) {
            minDiff = diff;
            closestPair = pair;
        }
    }
    return closestPair;
}

void printMatrix(const vector<vector<int>>& X, int n, ofstream& logfile) {
    int width = 6;  // spacing for alignment

    // Print column indices
    logfile << setw(width) << " ";
    for (int j = 0; j <= n; ++j)
        logfile << setw(width) << j;
    logfile << "\n";

    // Separator line
    logfile << setw(width) << " ";
    for (int j = 0; j <= n; ++j)
        logfile << string(width, '-');
    logfile << "\n";

    // Print rows
    for (int i = 0; i <= n; ++i) {
        logfile << setw(width) << i;  // row index
        for (int j = 0; j <= n; ++j) {
            if (X[i][j] > 0)
                logfile << setw(width) << "+" + to_string(X[i][j]);
            else
                logfile << setw(width) << X[i][j];
        }
        logfile << "\n";
    }
}

void printMatrixVisual(const vector<vector<int>>& X, int n, ofstream& logfile) {
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (X[i][j] != 0)
                logfile << "A";
            else
                logfile << "0";
        }
        logfile << "\n";
    }
}

// local alignment to decide rearrangement
string local_alignment(const string& seq, ofstream& logfile, int top_x = 1) {
    int n = seq.size();
    string rev_seq = reverse_complement(seq);

    const int NEG_INF = -1000000000;
    int match = 3, substitution = -5, gap_open = -6, gap_extend = -5;

    vector<vector<int>> dp(n + 1, vector<int>(n + 1, 0));
    vector<vector<int>> gap_in_seq1(n + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<int>> gap_in_seq2(n + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<pair<int, int>>> traceback(n + 1, vector<pair<int, int>>(n + 1, {-1, -1}));

    // Initialize first row and column for local alignment
    for (int i = 0; i <= n; ++i) {
        dp[i][0] = 0;
        gap_in_seq1[i][0] = NEG_INF;
        gap_in_seq2[i][0] = NEG_INF;
        traceback[i][0] = {-1, -1};
    }
    for (int j = 0; j <= n; ++j) {
        dp[0][j] = 0;
        gap_in_seq1[0][j] = NEG_INF;
        gap_in_seq2[0][j] = NEG_INF;
        traceback[0][j] = {-1, -1};
    }

    // DP fill
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            int scoreDiag = dp[i - 1][j - 1] + (seq[i - 1] == rev_seq[j - 1] ? match : substitution);
            int scoreLeft = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
            int scoreUp   = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

            dp[i][j] = max({0, scoreDiag, scoreLeft, scoreUp});
            gap_in_seq1[i][j] = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
            gap_in_seq2[i][j] = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

            if (dp[i][j] == 0)
                traceback[i][j] = {-1, -1};
            else if (dp[i][j] == scoreDiag)
                traceback[i][j] = {i - 1, j - 1};
            else if (dp[i][j] == scoreLeft)
                traceback[i][j] = {i, j - 1};
            else if (dp[i][j] == scoreUp)
                traceback[i][j] = {i - 1, j};
        }
    }

    // Collect all local maxima
    struct Cell { int i, j, score; };
    vector<Cell> candidates;

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (dp[i][j] <= 0) continue;
            bool is_max = true;
            for (int di = -1; di <= 1 && is_max; ++di) {
                for (int dj = -1; dj <= 1 && is_max; ++dj) {
                    if (di == 0 && dj == 0) continue;
                    int ni = i + di, nj = j + dj;
                    if (ni >= 1 && ni <= n && nj >= 1 && nj <= n) {
                        if (dp[i][j] < dp[ni][nj]) is_max = false;
                    }
                }
            }
            if (is_max) candidates.push_back({i, j, dp[i][j]});
        }
    }

    // Sort by descending score
    sort(candidates.begin(), candidates.end(), [](const Cell &a, const Cell &b) {
        return a.score > b.score;
    });

    vector<Cell> top_cells;
    vector<vector<int>> path_matrix(n + 1, vector<int>(n + 1, 0)); // merged matrix

    struct PathInfo {
        int start_i, start_j;
        int end_i, end_j;
        int score;
        pair<int,int> closest_pair; // store for rearrangement check
    };
    vector<PathInfo> path_infos;

    // Pick top-x distinct paths
    for (auto &cell : candidates) {
        if ((int)top_cells.size() >= top_x) break;

        int ti = cell.i, tj = cell.j;
        bool overlap = false;

        // Check for overlap
        while (ti > 0 && tj > 0 && dp[ti][tj] > 0) {
            if (path_matrix[ti][tj] > 0) { overlap = true; break; }
            auto [prev_i, prev_j] = traceback[ti][tj];
            if (prev_i == -1 || prev_j == -1) break;
            ti = prev_i; tj = prev_j;
        }

        if (!overlap) {
            // Trace again to mark and find start
            ti = cell.i; tj = cell.j;
            int start_i = ti, start_j = tj;
            while (ti > 0 && tj > 0 && dp[ti][tj] > 0) {
                path_matrix[ti][tj] = cell.score;
                auto [prev_i, prev_j] = traceback[ti][tj];
                if (prev_i == -1 || prev_j == -1 || dp[prev_i][prev_j] == 0) break;
                ti = prev_i; tj = prev_j;
            }
            start_i = ti; start_j = tj;

            // Temporary placeholder for pair
            path_infos.push_back({start_i, start_j, cell.i, cell.j, cell.score, {-1,-1}});
            top_cells.push_back(cell);
        }
    }

    logfile << "\n===== Local Alignment (Top " << top_x << ") paths =====\n";

    // Print alignments and store closest pairs
    for (auto &info : path_infos) {
        int i = info.end_i, j = info.end_j;
        int path_score = info.score;
        string aligned_seq1 = "", aligned_seq2 = "";

        vector<pair<int,int>> alignment_path;

        while (i > 0 && j > 0 && dp[i][j] > 0) {
            alignment_path.push_back({i - 1, j - 1});
            auto [prev_i, prev_j] = traceback[i][j];
            if (prev_i == -1 || prev_j == -1) break;

            if (prev_i == i - 1 && prev_j == j - 1) {
                aligned_seq1 = seq[i - 1] + aligned_seq1;
                aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
            } else if (prev_i == i && prev_j == j - 1) {
                aligned_seq1 = "-" + aligned_seq1;
                aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
            } else if (prev_i == i - 1 && prev_j == j) {
                aligned_seq1 = seq[i - 1] + aligned_seq1;
                aligned_seq2 = "-" + aligned_seq2;
            }

            if(dp[prev_i][prev_j] == 0) break;
            i = prev_i; j = prev_j;
        }

        // Reverse to print from start→end
        reverse(alignment_path.begin(), alignment_path.end());

        logfile << "Alignment score: " << path_score << "\n";
        logfile << "Start cell: (" << info.start_i - 1 << ", " << info.start_j - 1 << ")"
                << "  End cell: (" << info.end_i - 1 << ", " << info.end_j - 1 << ")\n";

        logfile << "Alignment path:\n";
        for (size_t k = 0; k < alignment_path.size(); ++k) {
            logfile << "(" << alignment_path[k].first << "," << alignment_path[k].second << ")";
            if (k + 1 < alignment_path.size()) logfile << "->";
        }
        logfile << "\n";

        vector<pair<int, int>> converted_alignment_path;
        for (const auto& [i, j] : alignment_path) {
            converted_alignment_path.push_back({i, n - (j) - 1});
        }

        logfile << "Converted alignment path:\n";
        for (size_t k = 0; k < converted_alignment_path.size(); ++k) {
            logfile << "(" << converted_alignment_path[k].first << "," << converted_alignment_path[k].second << ")";
            if (k + 1 < converted_alignment_path.size()) logfile << "->";
        }
        logfile << "\n";

        auto result = findClosestPair(converted_alignment_path);
        logfile << "Closest pair: (" << result.first << ", " << result.second << ")\n";
        info.closest_pair = result;  // store for rearrangement check

        logfile << aligned_seq1 << "\n" << aligned_seq2 << "\n\n";
    }

    // // Print combined matrix
    // logfile << "Top-" << top_x << " paths matrix:\n";
    // for (int i = 1; i <= n; ++i) {
    //     for (int j = 1; j <= n; ++j) {
    //         logfile << path_matrix[i][j] << "\t";
    //     }
    //     logfile << "\n";
    // }

    logfile << "====================================\n\n";

    // ===============================
    // 🔹 Rearrangement decision module
    // ===============================
    int gap_threshold = 2;
    logfile << "=== Rearrangement Decision Module ===\n";
    logfile << "Gap threshold = " << gap_threshold << "\n";

    for (size_t idx = 0; idx < path_infos.size(); ++idx) {
        auto [x, y] = path_infos[idx].closest_pair;
        int gap = y - x;
        logfile << "Path #" << idx + 1 << ": Closest pair (" << x << ", " << y
                << "), gap = " << gap << "\n";

        if (gap >= 0 && gap < gap_threshold) {
            logfile << "✅ Condition met (gap positive and < threshold)\n";
            logfile << "Rearranging sequence: right part [" << y << "..." << n-1
                    << "] + left part [0..." << y-1 << "]\n";

            string rearranged = seq.substr(y) + seq.substr(0, y);
            logfile << "Rearranged sequence: " << rearranged << "\n";
            logfile << "Stopping further checks (rearranged successfully).\n";
            logfile << "====================================\n\n";
            return rearranged;
        } else {
            if (gap <= 0)
                logfile << "❌ Condition not met (gap not positive)\n";
            else
                logfile << "❌ Condition not met (gap >= threshold)\n";
        }
    }

    logfile << "No rearrangement condition met. Returning original sequence.\n";
    logfile << "====================================\n\n";

    return seq;
}

// Perform modified global alignment between rearranged seq and its reverse complement
tuple<string, string, string, string> global_alignment_custom(const string& seq, ofstream& logfile) {
    int n = seq.size();
    string rev_seq = reverse_complement(seq);

    const int NEG_INF = -1000000000;
    int match = 3, substitution = -5, gap_open = -6, gap_extend = -5;

    vector<vector<int>> dp(n + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<int>> gap_in_seq1(n + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<int>> gap_in_seq2(n + 1, vector<int>(n + 1, NEG_INF));
    vector<vector<pair<int, int>>> traceback(n + 1, vector<pair<int, int>>(n + 1, {-1, -1}));

    // Initialize
    dp[0][0] = 0;
    for (int i = 1; i <= n; ++i) {
        dp[i][0] = gap_open + (i - 1) * gap_extend;
        traceback[i][0] = {i - 1, 0};
    }
    for (int j = 1; j <= n; ++j) {
        dp[0][j] = gap_open + (j - 1) * gap_extend;
        traceback[0][j] = {0, j - 1};
    }

    // Fill DP matrix
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            int scoreDiag = dp[i - 1][j - 1] + (seq[i - 1] == rev_seq[j - 1] ? match : substitution);
            int scoreLeft = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
            int scoreUp = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

            dp[i][j] = max({scoreDiag, scoreLeft, scoreUp});
            gap_in_seq1[i][j] = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
            gap_in_seq2[i][j] = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

            if (dp[i][j] == scoreDiag)
                traceback[i][j] = {i - 1, j - 1};
            else if (dp[i][j] == scoreLeft)
                traceback[i][j] = {i, j - 1};
            else
                traceback[i][j] = {i - 1, j};
        }
    }

    // Find best cell in region where i <= n - j
    int best_i = -1, best_j = -1, best_score = NEG_INF;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (i <= n - j && dp[i][j] > best_score) {
                best_score = dp[i][j];
                best_i = i;
                best_j = j;
            }
        }
    }

    // Traceback
    int i = best_i, j = best_j;
    string aligned_seq1 = "", aligned_seq2 = "";
    vector<pair<int, int>> alignment_path;

    while (i > 0 && j > 0) {
        alignment_path.push_back({i - 1, j - 1});
        auto [pi, pj] = traceback[i][j];
        if (pi == -1 || pj == -1) break;

        if (pi == i - 1 && pj == j - 1) {
            aligned_seq1 = seq[i - 1] + aligned_seq1;
            aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
        } else if (pi == i && pj == j - 1) {
            aligned_seq1 = "-" + aligned_seq1;
            aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
        } else if (pi == i - 1 && pj == j) {
            aligned_seq1 = seq[i - 1] + aligned_seq1;
            aligned_seq2 = "-" + aligned_seq2;
        }

        if (pi == i - 1 && pj == j - 1 && (i == 1 || j == 1)) break;
        i = pi;
        j = pj;
    }

    reverse(alignment_path.begin(), alignment_path.end());

    // Print all info like before
    logfile << "===== Global Alignment (custom) =====\n";
    logfile << "Max score cell (restricted area): (" << best_i - 1 << ", " << best_j - 1 << ") Score: " << best_score << "\n";
    logfile << "Start cell (traceback end): (" << i  - 1 << ", " << j - 1 << ")\n";
    logfile << "Normalized max score: " << (double) best_score / n << "\n";

    logfile << "Alignment path:\n";
    for (size_t k = 0; k < alignment_path.size(); ++k) {
        logfile << "(" << alignment_path[k].first << "," << alignment_path[k].second << ")";
        if (k + 1 < alignment_path.size()) logfile << "->";
    }
    logfile << "\n";

    logfile << aligned_seq1 << "\n" << aligned_seq2 << "\n";

    int first_cut = best_i - 1;
    int second_cut = n - (best_j - 1) - 1;

    logfile << "length of unit: " << n << "\n";
    logfile << "Cut 1: " << first_cut<< "\n";
    logfile << "Cut 2: " << second_cut << "\n";

    //Get A and partial A'
    // two options, A is longer of the two

    string right1 = seq.substr(first_cut + 1);       // from first_cut + 1 to end
    string left1 = seq.substr(0, first_cut + 1);    // from 0 to first_cut
    string left2 = seq.substr(0, second_cut + 1);      // from 0 to second_cut
    string right2 = seq.substr(second_cut + 1);        // from second_cut + 1 to end

    // Determine A1 and Abar1 based on which side is longer
    string A1, Abar1;
    if (right1.size() >= left1.size()) {
        A1 = right1;
        Abar1 = left1;
    } else {
        A1 = left1;
        Abar1 = right1;
    }

    // Determine A2 and Abar2 based on which side is longer
    string A2, Abar2;
    if (right2.size() >= left2.size()) {
        A2 = right2;
        Abar2 = left2;
    } else {
        A2 = left2;
        Abar2 = right2;
    }

    // Log everything
    logfile << "\n--- Cut 1 ---\n";
    logfile << "left1  (" << left1.size()  << "): " << left1  << "\n";
    logfile << "right1 (" << right1.size() << "): " << right1 << "\n";
    logfile << "A1     (" << A1.size()     << "): " << A1     << "\n";
    logfile << "Abar1  (" << Abar1.size()  << "): " << Abar1  << "\n";

    logfile << "\n--- Cut 2 ---\n";
    logfile << "left2  (" << left2.size()  << "): " << left2  << "\n";
    logfile << "right2 (" << right2.size() << "): " << right2 << "\n";
    logfile << "A2     (" << A2.size()     << "): " << A2     << "\n";
    logfile << "Abar2  (" << Abar2.size()  << "): " << Abar2  << "\n";

    double filter_threshold = 0.25;

    double ratio1 = (double) Abar1.size()/A1.size();
    double ratio2 = (double) Abar2.size()/A2.size();

    logfile << "ratio1 " << ratio1 << ", ratio2 " << ratio2 << endl;

    if(ratio1 < filter_threshold &&  ratio2 < filter_threshold)
    { 
        logfile << "\nAbar1 and Abar2 ratio " << ratio1 << " and " << ratio2 << " is smaller than " << filter_threshold << "\n";
        logfile << "====================================\n\n";
        return {"", "", "", ""};
    }

    if(ratio1 < filter_threshold)
    { 
        logfile << "\nAbar1 ratio " << ratio1 << " is smaller than " << filter_threshold << "\n";
        logfile << "====================================\n\n";
        return {"", "", A2, Abar2};
    }

    if(ratio2 < filter_threshold)
    { 
        logfile << "\nAbar2 ratio " << ratio2 << " is smaller than " << filter_threshold << "\n";
        logfile << "====================================\n\n";
        return {A1, Abar1, "", ""};
    }

    logfile << "====================================\n\n";

    return {A1, Abar1, A2, Abar2};
}


// // align extracted A to rev of pA' comp locally to see score
// double verify_AAbar_local(string A, string Abar, ofstream& logfile)
// {
//     string seq = A;
//     string rev_seq = reverse_complement(Abar);
//     int n = seq.size();
//     int m = rev_seq.size();

//     const int NEG_INF = -1000000000;
//     int match = 3, substitution = -5, gap_open = -6, gap_extend = -5;

//     vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
//     vector<vector<int>> gap_in_seq1(n + 1, vector<int>(m + 1, NEG_INF));
//     vector<vector<int>> gap_in_seq2(n + 1, vector<int>(m + 1, NEG_INF));
//     vector<vector<pair<int, int>>> traceback(n + 1, vector<pair<int, int>>(m + 1, {-1, -1}));

//     // Fill DP matrix (local alignment)
//     int max_score = 0;
//     int max_i = 0, max_j = 0;

//     for (int i = 1; i <= n; ++i) {
//         for (int j = 1; j <= m; ++j) {
//             int scoreDiag = dp[i - 1][j - 1] + (seq[i - 1] == rev_seq[j - 1] ? match : substitution);
//             int scoreLeft = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
//             int scoreUp = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

//             dp[i][j] = max({0, scoreDiag, scoreLeft, scoreUp});
//             gap_in_seq1[i][j] = max(dp[i][j - 1] + gap_open, gap_in_seq1[i][j - 1] + gap_extend);
//             gap_in_seq2[i][j] = max(dp[i - 1][j] + gap_open, gap_in_seq2[i - 1][j] + gap_extend);

//             if (dp[i][j] == 0) traceback[i][j] = {-1, -1};
//             else if (dp[i][j] == scoreDiag) traceback[i][j] = {i - 1, j - 1};
//             else if (dp[i][j] == scoreLeft) traceback[i][j] = {i, j - 1};
//             else if (dp[i][j] == scoreUp) traceback[i][j] = {i - 1, j};

//             if (dp[i][j] > max_score) {
//                 max_score = dp[i][j];
//                 max_i = i;
//                 max_j = j;
//             }
//         }
//     }

//     // Traceback from max cell
//     int i = max_i, j = max_j;
//     string aligned_seq1 = "", aligned_seq2 = "";
//     vector<pair<int,int>> alignment_path;

//     while (i > 0 && j > 0 && dp[i][j] > 0) {
//         alignment_path.push_back({i - 1, j - 1});
//         auto [pi, pj] = traceback[i][j];
//         if (pi == -1 || pj == -1 || dp[pi][pj] == 0) break;

//         if (pi == i - 1 && pj == j - 1) {
//             aligned_seq1 = seq[i - 1] + aligned_seq1;
//             aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
//         } else if (pi == i && pj == j - 1) {
//             aligned_seq1 = "-" + aligned_seq1;
//             aligned_seq2 = rev_seq[j - 1] + aligned_seq2;
//         } else if (pi == i - 1 && pj == j) {
//             aligned_seq1 = seq[i - 1] + aligned_seq1;
//             aligned_seq2 = "-" + aligned_seq2;
//         }

//         i = pi;
//         j = pj;
//     }

//     reverse(alignment_path.begin(), alignment_path.end());

//     // Print info similar to local_alignment_custom_local
//     logfile << "===== Local Alignment of A vs rev(Abar) =====\n";
//     logfile << "Max score cell: (" << max_i - 1 << ", " << max_j - 1 << ") Score: " << max_score << "\n";
//     logfile << "Start cell (traceback end): (" << i - 1 << ", " << j - 1 << ")\n";
//     logfile << "Normalized max score: " << (double)max_score / n << "\n";

//     logfile << "Alignment path:\n";
//     for (size_t k = 0; k < alignment_path.size(); ++k) {
//         logfile << "(" << alignment_path[k].first << "," << alignment_path[k].second << ")";
//         if (k + 1 < alignment_path.size()) logfile << "->";
//     }
//     logfile << "\n";

//     logfile << aligned_seq1 << "\n" << aligned_seq2 << "\n";
//     logfile << "Length of A: " << n << "\n";
//     logfile << "Length of pA': " << m << "\n";

//     logfile << "====================================\n\n";

//     return (double)max_score / n;
// }


int main(int argc, char* argv[]) {
    if (argc != 9) {
        cerr << "Usage: " << argv[0] << " <input_fasta> <A1_file> <Abar1_file> <A2_file> <Abar2_file> <AAbar_file> <rearranged_file> <log_file>" << endl;
        return 1;
    }

    string input_fasta = argv[1];
    string output_file_A1 = argv[2];
    string output_file_Abar1 = argv[3];
    string output_file_A2 = argv[4];
    string output_file_Abar2 = argv[5];
    string output_file_AAbar = argv[6];
    string rearranged_file = argv[7];
    string log_file = argv[8];

    ifstream infile(input_fasta);

    ofstream A1_file(output_file_A1);
    ofstream Abar1_file(output_file_Abar1);

    ofstream A2_file(output_file_A2);
    ofstream Abar2_file(output_file_Abar2);

    ofstream AAbar_file(output_file_AAbar);

    ofstream rearrangedfile(rearranged_file);
    ofstream logfile(log_file);

    if (!infile.is_open() || !A1_file.is_open() || !Abar1_file.is_open() || !A2_file.is_open() || !Abar2_file.is_open() || !AAbar_file.is_open() || !rearrangedfile.is_open() || !logfile.is_open()) {
        cerr << "Error opening one or more files." << endl;
        return 1;
    }

    string header, seq;

    while (getline(infile, header)) {

        getline(infile, seq);
        logfile << header << "\n";

        if(seq == "norepeat" || seq == "nocycle" || seq == "toolong")
        {
            logfile << "norepeat/nocycle/toolong\n\n";
            continue;
        }

        string rearranged_seq = local_alignment(seq, logfile, 3);

        // write rearranged sequence to file
        rearrangedfile << header << " length_rearranged=" << rearranged_seq.size() << "\n";
        rearrangedfile << rearranged_seq << "\n";

        auto [A1, Abar1, A2, Abar2]  = global_alignment_custom(rearranged_seq, logfile);

        if(A1 == "" && A2 == "")
        {
            logfile << "(empty A1 and A2)\n\n";
            continue;
        }

        // ----- Write A1 set only if A1 is non-empty -----
        if (!A1.empty()) {
            // Write to A1_file
            A1_file << header << " lengthA=" << A1.size() << "\n";
            A1_file << A1 << "\n";

            // Write to Abar1_file
            Abar1_file << header << " lengthAbar=" << Abar1.size() << "\n";
            Abar1_file << Abar1 << "\n";
        }

        // ----- Write A2 set only if A2 is non-empty -----
        if (!A2.empty()) {
            // Write to A2_file
            A2_file << header << " lengthA=" << A2.size() << "\n";
            A2_file << A2 << "\n";

            // Write to Abar2_file
            Abar2_file << header << " lengthAbar=" << Abar2.size() << "\n";
            Abar2_file << Abar2 << "\n";
        }

        // Write to AAbar_file
        AAbar_file << header << " lengthAAbar=" << seq.size() << "\n";
        AAbar_file << seq << "\n";
    }

    infile.close();
    A1_file.close();
    Abar1_file.close();
    A2_file.close();
    Abar2_file.close();
    AAbar_file.close();
    rearrangedfile.close();
    logfile.close();

    return 0;
}
