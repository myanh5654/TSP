#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm> 
#include <numeric>
#include <iomanip>
#include <fstream> 
#include <cstdlib> 
#include <ctime> 
#include <chrono> 

using namespace std;

// ============================================================
// 1. KHAI BÁO BIẾN TOÀN CỤC & CẤU TRÚC DỮ LIỆU
// ============================================================

int N = 0; 
vector<vector<double>> distMatrix; 
vector<vector<int>> tabuMatrix; 

struct City {
    int id;
    double x, y;
};
vector<City> cities; 

// ============================================================
// 2. CÁC HÀM HỖ TRỢ
// ============================================================

double euclidean_distance(const City& c1, const City& c2) {
    return sqrt(pow(c1.x - c2.x, 2) + pow(c1.y - c2.y, 2));
}

void precomputeDistances() {
    distMatrix.assign(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            distMatrix[i][j] = euclidean_distance(cities[i], cities[j]);
        }
    }
}

void readInput(const string& filename) {
    ifstream inFile(filename); 
    if (!inFile) { cerr << "LỖI: Không thể mở tệp " << filename << endl; exit(1); }

    string line;
    string edgeType = "EXPLICIT"; 
    bool inSection = false;
    cities.clear();
    vector<double> weights;

    while (getline(inFile, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n")); 
        if (line.empty()) continue;

        if (line.find("DIMENSION") != string::npos) {
            size_t num_start = line.find_first_of("0123456789");
            if (num_start != string::npos) N = stoi(line.substr(num_start));
        } 
        else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
            if (line.find("EUC_2D") != string::npos) edgeType = "EUC_2D";
            else if (line.find("EXPLICIT") != string::npos) edgeType = "EXPLICIT";
        }
        else if (line.find("NODE_COORD_SECTION") != string::npos) { inSection = true; continue; }
        else if (line.find("EDGE_WEIGHT_SECTION") != string::npos) { inSection = true; continue; }
        else if (line.find("EOF") != string::npos) break;

        if (inSection) {
            stringstream ss(line);
            if (edgeType == "EUC_2D") {
                int id; double x, y;
                if (ss >> id >> x >> y) cities.push_back({id - 1, x, y});
            } else {
                double val; while (ss >> val) weights.push_back(val);
            }
        }
    }
    inFile.close();

    if (N == 0) { cerr << "Lỗi DIMENSION!" << endl; exit(1); }

    if (edgeType == "EUC_2D") {
        if (cities.size() < N) N = cities.size();
        precomputeDistances(); 
    } else {
        distMatrix.assign(N, vector<double>(N, 0.0));
        int k = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) { 
                if (k < weights.size()) distMatrix[i][j] = distMatrix[j][i] = weights[k++];
            }
        }
    }
}

double calculateTourCost(const vector<int>& tour) {
    double cost = 0.0;
    for (int i = 0; i < N; ++i) {
        cost += distMatrix[tour[i]][tour[(i + 1) % N]];
    }
    return cost;
}

vector<int> createNearestNeighborTour() {
    vector<int> tour; 
    vector<bool> visited(N, false); 
    int cur = 0; 
    tour.push_back(cur); 
    visited[cur] = true;

    for (int i = 0; i < N - 1; ++i) { 
        double min_d = 1e9; 
        int next = -1;
        for (int j = 0; j < N; ++j) { 
            if (!visited[j] && distMatrix[cur][j] < min_d) { 
                min_d = distMatrix[cur][j]; 
                next = j; 
            }
        }
        if (next == -1) break; 
        cur = next; 
        tour.push_back(cur); 
        visited[cur] = true; 
    }
    return tour;
}

void diversifyTour(vector<int>& tour, int strength) {
    for (int i = 0; i < strength; ++i) { 
        int p1 = rand() % N; 
        int p2 = rand() % N; 
        if (p1 != p2) swap(tour[p1], tour[p2]); 
    }
}

// ============================================================
// 3. THUẬT TOÁN TABU SEARCH (FIXED LOGIC)
// ============================================================

string runTabuSearch(const string& instanceName, int denom, int maxIterations) {
    auto start_time = chrono::high_resolution_clock::now();

    const int TABU_TENURE = (N > 0 && denom > 0) ? (int)round((double)N / denom) : 10; 
    const int MAX_ITERATIONS = maxIterations; 
    const int DIVERSIFY_THRESHOLD = 150; 
    
    // --- Tạo file Log ---
    string logFilename = instanceName;
    size_t lastdot = logFilename.find_last_of(".");
    if (lastdot != string::npos) logFilename.erase(lastdot);
    logFilename += "_" + to_string(denom) + "_log.txt";
    ofstream logFile(logFilename);
    
    logFile << "Current_Iter,Solution_Cost,Best_Cost,Next_Best_Found,Move_Type,Move_Desc,Tabu_List" << endl;
    
    vector<int> currentTour = createNearestNeighborTour(); 
    double currentCost = calculateTourCost(currentTour);
    if (currentCost < 0) currentCost = abs(currentCost); 

    vector<int> bestTour = currentTour;
    double bestCost = currentCost;
    int bestIter = 0; 
    
    tabuMatrix.assign(N, vector<int>(N, 0)); 
    int noImprovement = 0; 

    // --- Vòng lặp chính ---
    for (int k = 1; k <= MAX_ITERATIONS; ++k) {
        noImprovement++;
        int moveType = rand() % 2; 
        double bestDelta = numeric_limits<double>::max();
        int best_i = -1, best_j = -1; 
        bool moveFound = false;

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;

                double delta = 0;
                bool isTabu = false;
                
                // === MOVE 1-0 (INSERTION) ===
                if (moveType == 0) {
                    if (i == j || i == (j + 1) % N) continue; 
                    int C = currentTour[i];
                    int A = currentTour[(i - 1 + N) % N];
                    int B = currentTour[(i + 1) % N];
                    int X = currentTour[j];
                    int Y = currentTour[(j + 1) % N];

                    double rem = distMatrix[A][C] + distMatrix[C][B] + distMatrix[X][Y];
                    double add = distMatrix[A][B] + distMatrix[X][C] + distMatrix[C][Y];
                    delta = add - rem;

                    // Kiểm tra các cạnh MỚI (A-B, X-C, C-Y) có bị cấm không?
                    if (tabuMatrix[A][B] > k || tabuMatrix[X][C] > k || tabuMatrix[C][Y] > k) {
                        isTabu = true;
                    }
                }
                // === MOVE 2-0 (BLOCK INSERTION) ===
                else if (moveType == 1) {
                    if (i >= N - 1) continue; 
                    int i_next = i + 1;
                    if (j == i || j == i_next || j == (i - 1 + N) % N) continue;

                    int C = currentTour[i];
                    int D = currentTour[i_next];
                    int A = currentTour[(i - 1 + N) % N];
                    int B = currentTour[(i + 2) % N];
                    int X = currentTour[j];
                    int Y = currentTour[(j + 1) % N];

                    double rem = distMatrix[A][C] + distMatrix[D][B] + distMatrix[X][Y];
                    double add = distMatrix[A][B] + distMatrix[X][C] + distMatrix[D][Y];
                    delta = add - rem;

                    // Kiểm tra các cạnh MỚI (A-B, X-C, D-Y) có bị cấm không?
                    if (tabuMatrix[A][B] > k || tabuMatrix[X][C] > k || tabuMatrix[D][Y] > k) {
                        isTabu = true;
                    }
                }

                if (currentCost + delta > 0) { 
                    bool aspiration = (currentCost + delta < bestCost); 
                    if (delta < bestDelta) {
                        if (!isTabu || aspiration) {
                            bestDelta = delta;
                            best_i = i;
                            best_j = j;
                            moveFound = true;
                        }
                    }
                }
            }
        }

        if (moveFound) {
            string moveDescStr = "";
            string moveTypeStr = "";
            string nextBestStr = (currentCost + bestDelta < bestCost) ? "Yes" : "No";

            if (moveType == 0) { 
                moveTypeStr = "1-0 (Ins)";
                int id_C = currentTour[best_i];
                int id_X = currentTour[best_j]; 
                
                int A = currentTour[(best_i - 1 + N) % N];
                int B = currentTour[(best_i + 1) % N];
                int Y = currentTour[(best_j + 1) % N]; // Thành phố sau X

                // [SỬA LỖI] Cập nhật Tabu: Cấm các cạnh VỪA BỊ XÓA (A-C, C-B, X-Y)
                // Để ngăn quay lại trạng thái cũ
                tabuMatrix[A][id_C] = tabuMatrix[id_C][A] = k + TABU_TENURE;
                tabuMatrix[id_C][B] = tabuMatrix[B][id_C] = k + TABU_TENURE;
                tabuMatrix[id_X][Y] = tabuMatrix[Y][id_X] = k + TABU_TENURE;

                moveDescStr = to_string(id_C + 1) + "-" + to_string(id_X + 1);

                currentTour.erase(currentTour.begin() + best_i);
                auto it = find(currentTour.begin(), currentTour.end(), id_X); 
                int new_pos = distance(currentTour.begin(), it);
                currentTour.insert(currentTour.begin() + new_pos + 1, id_C);
            } 
            else if (moveType == 1) { 
                moveTypeStr = "2-0 (Blk)";
                int i = best_i; 
                int j = best_j;
                
                int id_C = currentTour[i];
                int id_D = currentTour[i+1];
                int id_X = currentTour[j]; 
                
                int A = currentTour[(i - 1 + N) % N];
                int B = currentTour[(i + 2) % N];
                int Y = currentTour[(j + 1) % N]; // Thành phố sau X

                // [SỬA LỖI] Cập nhật Tabu: Cấm các cạnh VỪA BỊ XÓA (A-C, D-B, X-Y)
                tabuMatrix[A][id_C] = tabuMatrix[id_C][A] = k + TABU_TENURE;
                tabuMatrix[id_D][B] = tabuMatrix[B][id_D] = k + TABU_TENURE;
                tabuMatrix[id_X][Y] = tabuMatrix[Y][id_X] = k + TABU_TENURE;

                moveDescStr = "[" + to_string(id_C + 1) + "," + to_string(id_D + 1) + "]-" + to_string(id_X + 1);

                currentTour.erase(currentTour.begin() + i, currentTour.begin() + i + 2); 
                auto it = find(currentTour.begin(), currentTour.end(), id_X); 
                int pos_X = distance(currentTour.begin(), it);
                currentTour.insert(currentTour.begin() + pos_X + 1, id_D);
                currentTour.insert(currentTour.begin() + pos_X + 1, id_C);
            }

            currentCost += bestDelta;

            if (currentCost < bestCost) {
                bestCost = currentCost;
                bestTour = currentTour;
                bestIter = k; 
                noImprovement = 0; 
            }

            // --- IN LOG BẢNG VÀ DANH SÁCH TABU ---
            string tabuListStr = "";
            int countTabu = 0;
            // Duyệt ma trận đối xứng (chỉ cần tam giác trên)
            for(int r=0; r<N; ++r) { 
                for(int c=r+1; c<N; ++c) {
                    if(tabuMatrix[r][c] > k) { 
                        tabuListStr += "(" + to_string(r+1) + "-" + to_string(c+1) + ":" + to_string(tabuMatrix[r][c]) + ");";
                        countTabu++;
                    }
                    if(countTabu >= 8) break; 
                }
                if (countTabu >= 8) { tabuListStr += "..."; break; }
            }
            if (tabuListStr.empty()) tabuListStr = "EMPTY";

            logFile << fixed << setprecision(2) << k << "," << currentCost << "," << bestCost << "," << nextBestStr << "," << moveTypeStr << "," << moveDescStr << "," << tabuListStr << endl;

            // In ra màn hình mỗi 100 vòng hoặc khi có Next Best
            if (k % 100 == 0 || k == 1 || nextBestStr == "Yes") { 
                cout << left 
                     << setw(8) << k 
                     << setw(15) << fixed << setprecision(2) << currentCost 
                     << setw(15) << fixed << setprecision(2) << bestCost    
                     << setw(10) << nextBestStr 
                     << setw(12) << moveTypeStr
                     << setw(20) << moveDescStr << endl;
                     
                cout << "   |_ Tabu: " << tabuListStr << endl; // In danh sách cạnh bị cấm
            }
        }

        if (noImprovement >= DIVERSIFY_THRESHOLD) {
            logFile << "DIVERSIFYING at iteration " << k << endl;
            diversifyTour(currentTour, N); 
            currentCost = calculateTourCost(currentTour); 
            tabuMatrix.assign(N, vector<int>(N, 0)); 
            noImprovement = 0;
        }
    }
    
    logFile.close();

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;
    double runtime = elapsed.count();
    
    stringstream ss;
    ss << fixed << setprecision(2) << denom << "," << bestCost << "," << bestIter << "," << fixed << setprecision(4) << runtime;
    return ss.str();
}

// ============================================================
// 4. HÀM MAIN
// ============================================================
int main(int argc, char* argv[]) {
    srand(unsigned(time(0))); 
    
    if (argc != 4) {
        cerr << "Usage: ./tsp-solver <file> <denom> <iter>" << endl;
        return 1;
    }

    string filename = argv[1]; 
    int denom = stoi(argv[2]);
    int maxIterations = stoi(argv[3]);

    readInput(filename); 

    if (N > 0) {
        cout << runTabuSearch(filename, denom, maxIterations) << endl; 
    } 
    return 0;
}