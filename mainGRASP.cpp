#include "Interface.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

int main(){
    File file("rc204");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    cout << "Total profit is --> " << totalProfit << '\n';
    for(int i = 0; i < 5; i++){
        TOP top(file.get_N(), file.get_V());
        auto start = high_resolution_clock::now();
        vector<int> firstSolution = top.grasp_solution_generator();
        int rcl_len = file.get_N() / file.get_V();
        GRASP grasp(file.get_N(), file.get_V(), rcl_len, 100,30, firstSolution);
        grasp.graspAlgorithm();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<seconds>(stop - start);
        cout << "**************************************************" << '\n';
        cout << "Time taken: "
             << duration.count() << " seconds" << endl;

        cout<<"Iteration "<<i<<" finished!"<<'\n';
    }

    return 0;
}
