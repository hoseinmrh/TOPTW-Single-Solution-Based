#include "Interface.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

int main() {
// Testing
    File file("c204");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    for(int iter = 0; iter< 8; iter++) {
        cout << "Total profit is --> " << totalProfit << '\n';
        TOP top(file.get_N(), file.get_V());
        vector<int> firstSolution = top.random_solution_generator();
        cout << "***************** First Solution *****************" << '\n';
        for (int i = 0; i < firstSolution.size(); i++) {
            cout << firstSolution[i] << ' ';
        }
        cout << '\n';
        cout << "**************************************************" << '\n';
        auto start = high_resolution_clock::now();
        SA sa(firstSolution, 0.3, 0.999, 30, 4000, file.get_N(), file.get_V());
        float finalProfit = sa.saAlgorithm();
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<seconds>(stop - start);
        cout << "**************************************************" << '\n';
        cout << "Time taken: "
             << duration.count() << " seconds" << endl;

        cout << "Accuracy is --> " << (finalProfit / totalProfit) << '\n';

        cout<<"Iteration "<<iter+1<<" finished"<<'\n';
        cout<<'\n';
    }
    return 0;
}
