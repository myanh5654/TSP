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
    // Đoạn code đọc file giữ nguyên từ file gốc
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
// 3. THUẬT TOÁN TABU SEARCH (FULL LOG VERSION)
// ============================================================

// Hàm này trả về một chuỗi chứa kết quả tổng kết (Instance, Cost, Time, Denom)
string runTabuSearch(const string& instanceName, int denom, int maxIterations) {
    auto start_time = chrono::high_resolution_clock::now();

    // Tính toán Tabu Tenure dựa trên N và Denominator (Hệ số Tabu)
    const int TABU_TENURE = (N > 0 && denom > 0) ? (int)round((double)N / denom) : 10; 
    const int MAX_ITERATIONS = maxIterations; 
    const int DIVERSIFY_THRESHOLD = 150; 
    
    // Tạo tên file log chi tiết: [instance]_[denom]_log.txt
    string logFilename = instanceName;
    size_t lastdot = logFilename.find_last_of(".");
    if (lastdot != string::npos) logFilename.erase(lastdot);
    logFilename += "_" + to_string(denom) + "_log.txt";
    ofstream logFile(logFilename);
    
    // Ghi tiêu đề log chi tiết
    logFile << "Current_Iter,Solution_Cost,Best_Cost,Next_Best_Found,Move_Type,Move_Desc,Tabu_List" << endl;
    
    vector<int> currentTour = createNearestNeighborTour(); 
    double currentCost = calculateTourCost(currentTour);
    if (currentCost < 0) currentCost = abs(currentCost); 

    vector<int> bestTour = currentTour;
    double bestCost = currentCost;
    int bestIter = 0; 
    
    tabuMatrix.assign(N, vector<int>(N, 0)); 
    int noImprovement = 0; 

    // Vòng lặp chính
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
                
                // Move 1-0 (Insert)
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
                    // Kiểm tra Tabu: Nếu bất kỳ cạnh nào liên quan đến việc di chuyển C bị cấm
                    // Ta chỉ cấm các cạnh (A, C), (C, B), (X, C), (C, Y)
                    if (tabuMatrix[A][C] > k || tabuMatrix[C][B] > k || tabuMatrix[X][C] > k || tabuMatrix[C][Y] > k) isTabu = true;
                }
                // Move 2-0 (Block Insert)
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
                    // Kiểm tra Tabu: Chỉ cấm các cạnh liên quan đến khối (C, D)
                    if (tabuMatrix[A][C] > k || tabuMatrix[D][B] > k || tabuMatrix[X][C] > k || tabuMatrix[D][Y] > k) isTabu = true;
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

            // Thực hiện di chuyển và cập nhật Tabu
            if (moveType == 0) { 
                moveTypeStr = "1-0 (Ins)";
                int id_C = currentTour[best_i];
                int X_pos = best_j; 
                
                // Loại bỏ C
                currentTour.erase(currentTour.begin() + best_i);
                
                // Vị trí chèn C (sau X)
                if (best_i < X_pos) X_pos--; // Bù trừ do erase
                
                currentTour.insert(currentTour.begin() + X_pos + 1, id_C);

                // Cập nhật Tabu (C bị cấm di chuyển)
                tabuMatrix[id_C][id_C] = k + TABU_TENURE; // Đơn giản hóa: cấm node C di chuyển
                moveDescStr = to_string(id_C + 1) + " to pos after " + to_string(currentTour[X_pos] + 1);
            } 
            else if (moveType == 1) { 
                moveTypeStr = "2-0 (Blk)";
                int i = best_i; 
                int i_next = i + 1;
                int j = best_j; // Chỉ mục của X
                int id_C = currentTour[i];
                int id_D = currentTour[i_next];
                int id_X = currentTour[j]; 

                // Lấy khối (C, D)
                vector<int> block = {currentTour[i], currentTour[i_next]};
                
                // Loại bỏ khối
                currentTour.erase(currentTour.begin() + i, currentTour.begin() + i + 2); 
                
                // Vị trí chèn khối (sau X)
                if (j > i) j -= 2; // Bù trừ do erase 2 phần tử
                
                currentTour.insert(currentTour.begin() + j + 1, block.begin(), block.end());
                
                // Cập nhật Tabu (C và D bị cấm di chuyển)
                tabuMatrix[id_C][id_C] = k + TABU_TENURE; 
                tabuMatrix[id_D][id_D] = k + TABU_TENURE; 
                
                moveDescStr = "[" + to_string(id_C + 1) + "," + to_string(id_D + 1) + "] to pos after " + to_string(id_X + 1);
            }

            currentCost += bestDelta;

            if (currentCost < bestCost) {
                bestCost = currentCost;
                bestTour = currentTour;
                bestIter = k; 
                noImprovement = 0; 
            }

            // --- GHI LOG CHI TIẾT ---
            // 1. Chuẩn bị Tabu List string
            string tabuListStr = "";
            int countTabu = 0;
            for(int r=0; r<N; ++r) { 
                if(tabuMatrix[r][r] > k) { // Chỉ cấm node, nên chỉ check chéo
                    tabuListStr += "(" + to_string(r+1) + ":" + to_string(tabuMatrix[r][r]) + ");";
                    countTabu++;
                }
                if (countTabu > 10) { tabuListStr += "..."; break; }
            }
            if (tabuListStr.empty()) tabuListStr = "EMPTY";

            // 2. Ghi dòng log
            logFile << fixed << setprecision(2) 
                    << k << ","
                    << currentCost << ","
                    << bestCost << ","
                    << nextBestStr << ","
                    << moveTypeStr << ","
                    << moveDescStr << ","
                    << tabuListStr << endl;

            // --- IN BẢNG TRÊN CONSOLE (để tiện debug/xem nhanh) ---
            if (k % 100 == 0 || k == 1) { // In mỗi 100 lần lặp hoặc lần đầu
                cout << left 
                     << setw(8) << k 
                     << setw(15) << fixed << setprecision(2) << currentCost 
                     << setw(15) << fixed << setprecision(2) << bestCost    
                     << setw(10) << nextBestStr 
                     << setw(12) << moveTypeStr
                     << setw(20) << moveDescStr << endl;
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
    
    // --- TRẢ VỀ CHUỖI KẾT QUẢ TỔNG KẾT (CSV) ---
    // Định dạng: Denom,BestCost,BestIter,Runtime
    stringstream ss;
    ss << fixed << setprecision(2) << denom << ","
       << bestCost << ","
       << bestIter << ","
       << fixed << setprecision(4) << runtime;
    return ss.str();
}

// ============================================================
// 4. HÀM MAIN
// ============================================================
int main(int argc, char* argv[]) {
    srand(unsigned(time(0))); 
    
    // YÊU CẦU: ./tsp-solver [instance_file] [denom] [max_iter]
    if (argc != 4) {
        cerr << "Cách sử dụng: ./tsp-solver <instance_file> <tabu_denominator> <max_iterations>" << endl;
        cerr << "Ví dụ: ./tsp-solver eil51.tsp 10 2000" << endl;
        return 1;
    }

    string filename = argv[1]; 
    int denom = stoi(argv[2]);
    int maxIterations = stoi(argv[3]);

    readInput(filename); 

    if (N > 0) {
        // In ra kết quả tổng kết cuối cùng (dạng CSV)
        cout << runTabuSearch(filename, denom, maxIterations) << endl; 
    } 
    return 0;
}