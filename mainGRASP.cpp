#include "Interface.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

int main(){
    File file("pr01");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    cout << "Total profit is --> " << totalProfit << '\n';
    TOP top(file.get_N(), file.get_V());
    vector<int> firstSolution = top.grasp_solution_generator();
//    cout << "***************** First Solution *****************" << '\n';
//    for (int i = 0; i < firstSolution.size(); i++) {
//        cout << firstSolution[i] << ' ';
//    }
    GRASP grasp(10);
    grasp.testRCL();
    cout << '\n';
}
