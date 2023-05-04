//
// Created by Hosein on 4/25/2023.
//

#include "Interface.h"
using namespace std;

int main(){
    File file("3-pr13");
    file.read_file();
    float totalProfit = 0;
    totalProfit = calculate_profit();
    vector<int> checkingVector = {0 ,142, 64, 48, 38, 7 ,93, 108, 70, 79, 112, 133, 0, -1, 0 ,82 ,123, 107, 87, 26, 139, 10 ,18, 25, 43, 84 ,23 ,130, 86, 16, 104, 58, 0 ,-1,0 ,138 ,128 ,33 ,37 ,24 ,19 ,78 ,5 ,22 ,4 ,117 ,6 ,61, 34, 141, 74 ,0,-1, 0, 110, 14 ,80, 52, 119, 27, 45 ,91 ,39 ,17 ,35 ,28 ,81, 44 ,31 ,144 ,0};
    TOP top(0);
    top.calculate_solution_final(checkingVector,1);
    cout<<"Finish"<<'\n';
    return 0;
}
