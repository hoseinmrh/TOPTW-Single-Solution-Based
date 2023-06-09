#include "Interface.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

int main(){
    File file("3-pr13");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    cout << "Total profit is --> " << totalProfit << '\n';
    for(int i = 0; i < 2; i++){
    TOP top(file.get_N(), file.get_V());
    auto start = high_resolution_clock::now();
    vector<int> firstSolution = top.random_solution_generator();
    TabuSearch tabuSearch (10, firstSolution, 5000);
    tabuSearch.tabuSearchAlgorithmTime(180);

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);
    cout << "**************************************************" << '\n';
    cout << "Time taken: "
         << duration.count() << " seconds" << endl;

        cout<<"Iteration "<<i<<" finished!"<<'\n';

        cout<<'\n';
    }

    return 0;
}