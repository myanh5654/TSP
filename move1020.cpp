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

using namespace std;

// ============================================================
// 1. KHAI BÁO BIẾN TOÀN CỤC & CẤU TRÚC DỮ LIỆU
// ============================================================

// N: Số lượng thành phố (được đọc từ file input)
int N = 0; 

// distMatrix: Ma trận lưu khoảng cách giữa mọi cặp thành phố
// distMatrix[i][j] là khoảng cách từ thành phố i đến j
vector<vector<double>> distMatrix; 

// tabuMatrix: Ma trận lưu trạng thái cấm (Tabu List)
// tabuMatrix[i][j] = k nghĩa là cạnh nối i và j bị cấm sử dụng cho đến vòng lặp thứ k
vector<vector<int>> tabuMatrix; 

// Cấu trúc lưu tọa độ của một thành phố
struct City {
    int id;      // ID của thành phố
    double x, y; // Tọa độ không gian
};
vector<City> cities; // Danh sách lưu tất cả các thành phố

// ============================================================
// 2. CÁC HÀM HỖ TRỢ TÍNH TOÁN & ĐỌC FILE
// ============================================================

/**
 * Tính khoảng cách Euclid giữa 2 thành phố dựa trên tọa độ
 */
double euclidean_distance(const City& c1, const City& c2) {
    return sqrt(pow(c1.x - c2.x, 2) + pow(c1.y - c2.y, 2));
}

/**
 * Tính toán trước khoảng cách giữa mọi cặp thành phố và lưu vào distMatrix
 * Giúp tra cứu nhanh (O(1)) trong quá trình chạy thuật toán
 */
void precomputeDistances() {
    // Khởi tạo ma trận kích thước N x N
    distMatrix.assign(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Tính khoảng cách và lưu vào ma trận
            distMatrix[i][j] = euclidean_distance(cities[i], cities[j]);
        }
    }
}

/**
 * Hàm đọc dữ liệu từ file .tsp
 * Hỗ trợ 2 định dạng: EUC_2D (Tọa độ) và EXPLICIT (Ma trận trọng số)
 */
void readInput(const string& filename) {
    ifstream inFile(filename); 
    if (!inFile) { cerr << "LỖI: Không thể mở tệp " << filename << endl; exit(1); }

    string line;
    string edgeType = "EXPLICIT"; // Mặc định kiểu dữ liệu
    bool inSection = false;       // Cờ đánh dấu đã vào phần dữ liệu chưa
    cities.clear();
    vector<double> weights;       // Dùng tạm để lưu trọng số nếu là file EXPLICIT

    // Đọc từng dòng của file
    while (getline(inFile, line)) {
        // Xóa khoảng trắng thừa đầu dòng
        line.erase(0, line.find_first_not_of(" \t\r\n")); 
        if (line.empty()) continue;

        // Tìm thông tin số lượng thành phố (DIMENSION)
        if (line.find("DIMENSION") != string::npos) {
            size_t num_start = line.find_first_of("0123456789");
            if (num_start != string::npos) N = stoi(line.substr(num_start));
        } 
        // Tìm kiểu dữ liệu (Tọa độ hay Ma trận)
        else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
            if (line.find("EUC_2D") != string::npos) edgeType = "EUC_2D";
            else if (line.find("EXPLICIT") != string::npos) edgeType = "EXPLICIT";
        }
        // Tìm điểm bắt đầu phần dữ liệu
        else if (line.find("NODE_COORD_SECTION") != string::npos) { inSection = true; continue; }
        else if (line.find("EDGE_WEIGHT_SECTION") != string::npos) { inSection = true; continue; }
        else if (line.find("EOF") != string::npos) break;

        // Đọc dữ liệu thực tế
        if (inSection) {
            stringstream ss(line);
            if (edgeType == "EUC_2D") {
                // Đọc ID, X, Y
                int id; double x, y;
                if (ss >> id >> x >> y) cities.push_back({id - 1, x, y});
            } else {
                // Đọc các con số trọng số liên tiếp
                double val; while (ss >> val) weights.push_back(val);
            }
        }
    }
    inFile.close();

    if (N == 0) { cerr << "Lỗi DIMENSION!" << endl; exit(1); }

    // Xử lý dữ liệu sau khi đọc xong
    if (edgeType == "EUC_2D") {
        if (cities.size() < N) N = cities.size();
        precomputeDistances(); // Tính khoảng cách từ tọa độ
    } else {
        // Chuyển mảng 1 chiều weights thành ma trận 2 chiều distMatrix
        distMatrix.assign(N, vector<double>(N, 0.0));
        int k = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j <= i; ++j) { // Thường là ma trận tam giác dưới
                if (k < weights.size()) distMatrix[i][j] = distMatrix[j][i] = weights[k++];
            }
        }
    }
}

/**
 * Tính tổng chi phí (độ dài) của một lộ trình (tour)
 */
double calculateTourCost(const vector<int>& tour) {
    double cost = 0.0;
    for (int i = 0; i < N; ++i) {
        // Cộng khoảng cách từ thành phố i đến thành phố tiếp theo (i+1)%N
        cost += distMatrix[tour[i]][tour[(i + 1) % N]];
    }
    return cost;
}

/**
 * Tạo lộ trình ban đầu bằng thuật toán Tham lam (Nearest Neighbor)
 * Giúp có một khởi đầu tốt hơn so với ngẫu nhiên
 */
vector<int> createNearestNeighborTour() {
    vector<int> tour; 
    vector<bool> visited(N, false); // Đánh dấu thành phố đã đi qua
    
    int cur = 0; // Bắt đầu từ thành phố 0
    tour.push_back(cur); 
    visited[cur] = true;

    for (int i = 0; i < N - 1; ++i) { 
        double min_d = 1e9; // Khởi tạo khoảng cách nhỏ nhất là vô cùng
        int next = -1;
        
        // Tìm thành phố gần nhất chưa đi qua
        for (int j = 0; j < N; ++j) { 
            if (!visited[j] && distMatrix[cur][j] < min_d) { 
                min_d = distMatrix[cur][j]; 
                next = j; 
            }
        }
        if (next == -1) break; 
        
        // Di chuyển đến thành phố tiếp theo
        cur = next; 
        tour.push_back(cur); 
        visited[cur] = true; 
    }
    return tour;
}

/**
 * Hàm Đa dạng hóa (Diversification)
 * Xáo trộn ngẫu nhiên lộ trình để thoát khỏi tối ưu cục bộ khi bị kẹt
 */
void diversifyTour(vector<int>& tour, int strength) {
    for (int i = 0; i < strength; ++i) { 
        int p1 = rand() % N; 
        int p2 = rand() % N; 
        if (p1 != p2) swap(tour[p1], tour[p2]); // Đổi chỗ ngẫu nhiên
    }
}

// ============================================================
// 3. THUẬT TOÁN TABU SEARCH (MOVE 1-0 & MOVE 2-0)
// ============================================================

void runTabuSearch() {
    // --- Cấu hình tham số ---
    const int MAX_ITERATIONS = 10000; // Số vòng lặp tối đa
    const int TABU_TENURE = N/5; // Thời gian cấm (số vòng lặp)
    const int DIVERSIFY_THRESHOLD = 150; // Ngưỡng kích hoạt đa dạng hóa (nếu bị kẹt)

    cout << "--- TABU SEARCH: MIXED MOVE 1-0 & 2-0 ---" << endl;
    
    // Tạo giải pháp ban đầu
    vector<int> currentTour = createNearestNeighborTour(); 
    double currentCost = calculateTourCost(currentTour);
    if (currentCost < 0) currentCost = abs(currentCost); // Đảm bảo chi phí dương

    // Lưu giữ giải pháp tốt nhất tìm được (Kỷ lục)
    vector<int> bestTour = currentTour;
    double bestCost = currentCost;
    
    // Khởi tạo Tabu Matrix (Ban đầu chưa cấm gì)
    tabuMatrix.assign(N, vector<int>(N, 0)); 
    int noImprovement = 0; // Đếm số vòng lặp chưa phá được kỷ lục

    cout << "[INIT] Cost: " << fixed << setprecision(2) << currentCost << endl;

    // --- Vòng lặp chính ---
    for (int k = 1; k <= MAX_ITERATIONS; ++k) {
        noImprovement++;

        // Random chọn loại Move cho vòng lặp này: 
        // 0: Move 1-0 (Chèn 1 thành phố)
        // 1: Move 2-0 (Chèn khối 2 thành phố)
        int moveType = rand() % 2; 

        double bestDelta = 1e9; // Delta tốt nhất trong các hàng xóm
        int best_i = -1, best_j = -1; // Lưu vị trí thực hiện move tốt nhất
        bool moveFound = false;

        // Duyệt qua tất cả các cặp vị trí (i, j) để tìm nước đi tốt nhất
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;

                double delta = 0;
                bool isTabu = false;
                
                // === TRƯỜNG HỢP 1: MOVE 1-0 (INSERTION) ===
                // Lấy thành phố tại i chèn vào sau j
                if (moveType == 0) {
                    if (i == j || i == (j + 1) % N) continue; // Bỏ qua nếu vị trí không đổi
                    
                    // Xác định các thành phố liên quan
                    int C = currentTour[i];                 // Thành phố bị di chuyển
                    int A = currentTour[(i - 1 + N) % N];   // Thành phố đứng trước C
                    int B = currentTour[(i + 1) % N];       // Thành phố đứng sau C
                    int X = currentTour[j];                 // Thành phố đích (chèn sau X)
                    int Y = currentTour[(j + 1) % N];       // Thành phố đứng sau X

                    // Tính toán sự thay đổi chi phí (Delta)
                    // Giảm: cạnh (A,C), (C,B), (X,Y)
                    // Tăng: cạnh (A,B), (X,C), (C,Y)
                    double rem = distMatrix[A][C] + distMatrix[C][B] + distMatrix[X][Y];
                    double add = distMatrix[A][B] + distMatrix[X][C] + distMatrix[C][Y];
                    delta = add - rem;

                    // Kiểm tra Tabu: Có cạnh mới nào đang bị cấm không?
                    if (tabuMatrix[A][B] > k || tabuMatrix[X][C] > k || tabuMatrix[C][Y] > k) isTabu = true;
                }

                // === TRƯỜNG HỢP 2: MOVE 2-0 (BLOCK INSERTION) ===
                // Lấy khối 2 thành phố [i, i+1] chèn vào sau j
                else if (moveType == 1) {
                    // Chỉ thực hiện nếu không phải cuối mảng (để đơn giản hóa logic index)
                    if (i >= N - 1) continue; 

                    int i_next = i + 1;
                    // Bỏ qua nếu vị trí đích j nằm trong khối hoặc ngay cạnh khối
                    if (j == i || j == i_next || j == (i - 1 + N) % N) continue;

                    int C = currentTour[i];         // Đầu khối
                    int D = currentTour[i_next];    // Cuối khối
                    int A = currentTour[(i - 1 + N) % N]; // Trước khối
                    int B = currentTour[(i + 2) % N];     // Sau khối
                    
                    int X = currentTour[j];         // Chèn sau X
                    int Y = currentTour[(j + 1) % N];

                    // Tính Delta cho khối
                    double rem = distMatrix[A][C] + distMatrix[D][B] + distMatrix[X][Y];
                    double add = distMatrix[A][B] + distMatrix[X][C] + distMatrix[D][Y];
                    delta = add - rem;

                    // Kiểm tra Tabu
                    if (tabuMatrix[A][B] > k || tabuMatrix[X][C] > k || tabuMatrix[D][Y] > k) isTabu = true;
                }

                // --- ĐÁNH GIÁ NƯỚC ĐI ---
                // Logic Tabu: Chỉ chấp nhận nếu KHÔNG bị cấm HOẶC Phá vỡ kỷ lục (Aspiration)
                if (currentCost + delta > 0) { // Đảm bảo không lỗi số âm
                    bool aspiration = (currentCost + delta < bestCost); // Tiêu chí khát vọng
                    
                    // Nếu tìm thấy delta tốt hơn delta tốt nhất hiện tại
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

        // --- THỰC HIỆN MOVE TỐT NHẤT TÌM ĐƯỢC ---
        if (moveFound) {
            int A, B; 
            
            // Thực hiện Move 1-0
            if (moveType == 0) { 
                int id_C = currentTour[best_i];
                int id_X = currentTour[best_j]; 
                
                // Cập nhật Tabu: Cấm quay lại trạng thái cũ
                A = currentTour[(best_i - 1 + N) % N];
                B = currentTour[(best_i + 1) % N];
                tabuMatrix[A][id_C] = tabuMatrix[id_C][A] = k + TABU_TENURE;
                tabuMatrix[id_C][B] = tabuMatrix[B][id_C] = k + TABU_TENURE;

                // Thao tác trên Vector: Xóa C và chèn vào sau X
                currentTour.erase(currentTour.begin() + best_i);
                auto it = find(currentTour.begin(), currentTour.end(), id_X); // Tìm lại vị trí X
                int new_pos = distance(currentTour.begin(), it);
                currentTour.insert(currentTour.begin() + new_pos + 1, id_C);
                
                cout << "Iter " << k << " [1-0]: Move " << id_C+1 << " after " << id_X+1 << " (Delta: " << bestDelta << ")" << endl;
            
            } 
            // Thực hiện Move 2-0
            else if (moveType == 1) { 
                int i = best_i; 
                int i_next = i + 1;
                int j = best_j;

                int id_C = currentTour[i];
                int id_D = currentTour[i_next];
                int id_X = currentTour[j]; 

                A = currentTour[(i - 1 + N) % N];
                // Cập nhật Tabu
                tabuMatrix[A][id_C] = tabuMatrix[id_C][A] = k + TABU_TENURE;

                // Thao tác trên Vector: Xóa khối [C,D] và chèn sau X
                currentTour.erase(currentTour.begin() + i, currentTour.begin() + i + 2); // Xóa 2 phần tử

                auto it = find(currentTour.begin(), currentTour.end(), id_X); // Tìm lại vị trí X
                int pos_X = distance(currentTour.begin(), it);
                
                // Chèn lại D rồi C (để thứ tự là C->D)
                currentTour.insert(currentTour.begin() + pos_X + 1, id_D);
                currentTour.insert(currentTour.begin() + pos_X + 1, id_C);

                cout << "Iter " << k << " [2-0]: Move Block [" << id_C+1 << "," << id_D+1 << "] after " << id_X+1 << " (Delta: " << bestDelta << ")" << endl;
            }

            // Cập nhật chi phí hiện tại
            currentCost += bestDelta;

            // Kiểm tra xem có phá kỷ lục không
            if (currentCost < bestCost) {
                bestCost = currentCost;
                bestTour = currentTour;
                noImprovement = 0; // Reset bộ đếm thất bại
                cout << "   *** NEW BEST! Cost: " << bestCost << " ***" << endl;
            }
        }

        // --- ĐA DẠNG HÓA (DIVERSIFICATION) ---
        // Nếu bị kẹt quá lâu (150 lần không cải thiện), xáo trộn để tìm hướng mới
        if (noImprovement >= DIVERSIFY_THRESHOLD) {
            cout << "--- STUCK! Diversifying... ---" << endl;
            diversifyTour(currentTour, N); // Xáo trộn mạnh
            currentCost = calculateTourCost(currentTour); // Tính lại chi phí
            tabuMatrix.assign(N, vector<int>(N, 0)); // Xóa sạch lệnh cấm cũ
            noImprovement = 0;
        }
    }

    // In kết quả cuối cùng
    cout << "\n--- FINAL RESULT ---" << endl;
    cout << "Final Best Cost: " << fixed << setprecision(2) << bestCost << endl;
    cout << "Final Best Tour: ";
    for(int i : bestTour) cout << i + 1 << " ";
    cout << endl;
}

// ============================================================
// 4. HÀM MAIN
// ============================================================
int main() {
    srand(unsigned(time(0))); // Khởi tạo bộ sinh số ngẫu nhiên
    string filename = "eil51.tsp";  // Tên file input

    readInput(filename); // Đọc dữ liệu
    
    if (N > 0) runTabuSearch(); // Chạy thuật toán
    
    return 0;
}