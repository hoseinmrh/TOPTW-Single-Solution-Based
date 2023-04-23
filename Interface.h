//
// Created by Hosein on 4/6/2023.
//

#ifndef HW_1_INTERFACE_H
#define HW_1_INTERFACE_H
#include <time.h>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>
#define MAX 20
using namespace std;

class File;
class Vertex;
class TOP;
class SA;
class GRASP;
class TabuItem;
static vector<Vertex> vertexVector;

void add_to_vertex_vector(Vertex v);

class Vertex{
protected:
    int i; //vertex number
    float x; //x coordinate
    float y; //y coordinate
    float d; //service duration
    float profit; //profit of the duration
    float opening_time;
    float closing_time;

    int nothing;

public:
    Vertex(int a){
        nothing = a;
    }
    Vertex(int i_v, float x_v, float y_v, float d_v, float profit_v, float openingTime_v, float closingTime_v){
        i = i_v;
        x = x_v;
        y = y_v;
        d = d_v;
        profit = profit_v;
        opening_time = openingTime_v;
        closing_time = closingTime_v;
    }

    int getI() const {
        return i;
    }

    float getX() const {
        return x;
    }

    float getY() const {
        return y;
    }

    float getD() const {
        return d;
    }

    float getProfit() const {
        return profit;
    }

    float getOpeningTime() const {
        return opening_time;
    }

    float getClosingTime() const {
        return closing_time;
    }

};

class File {
protected:
    int N, V;
    string fileName;
public:

    File(string file_name){
        fileName = file_name;
    }

    int get_N(){
        return N;
    }

    int get_V(){
        return V;
    }



    void set_N(int n){
        N = n;
    }

    void set_V(int v){
        V = v;
    }

    void read_file() {
        string firstPath = "C:\\HOSEIN\\jozve and tamrin\\8th Semester\\ADA\\HW 1\\Instances\\";
        string format = ".txt";
        string filePath = firstPath + fileName + format;
        ifstream new_file;
        // Open a file to perform a write operation using a file object.
        cout << filePath << '\n';
        new_file.open(filePath);

        if (!new_file) {
            cout << "Can't open file!" << '\n';
            exit(1);
        }
        if (new_file.is_open()) {
            string sa;
            int lineCount = 1;
            // Read data from the file object and put it into a string.
            while (getline(new_file, sa)) {
                if (lineCount == 2) {
                    lineCount++;
                    continue;
                }
                split_lines(sa, lineCount);
                lineCount++;
            }
            new_file.close();

        }
    }
    void split_lines(string line, int lineCount){
        float array[MAX];
        fill_n(array,MAX,-1);
        if(lineCount == 1){

            stringstream lineStream(line);
            int index = 0;
            while (lineStream.good() && index < MAX)
            {

                lineStream >> array[index];
                index ++;
            }
            set_V(array[1]);
            set_N(array[2]);

        }
        else{
            stringstream lineStream(line);
            int index = 0;
            int row = lineCount - 3;
            while (lineStream.good() && index < MAX)
            {


                lineStream >> array[index];
                index ++;
            }
            create_vertex(array);

        }
    }

    void create_vertex(float array[]){
        int i = array[0];
        float x = array[1];
        float y = array[2];
        float d = array[3];
        float profit = array[4];
        int c_index = not_minus1_index(array);
        float opening_time = array[c_index - 1];
        float closing_time = array[c_index];
        Vertex vertex(i,x,y,d,profit,opening_time,closing_time);
        add_to_vertex_vector(vertex);

    }

    int not_minus1_index(float array[]) {
        for (int i = MAX - 1; i > -1; i--) {
            if (array[i] != -1)
                return (i);

        }
    }


};

class TOP{
protected:
    int N;
    int V;
    float time;
    float profit;
    int nothing;

public:
    TOP(int a ){ //Temp Constructor
        nothing = a;
        time = 0;
        profit = 0;
    }

    TOP(int n, int v){
        N = n;
        V = v;
        time = 0;
        profit = 0;
    }

    float getProfit() const {
        return profit;
    }

    void addProfit(float profit){
        TOP::profit += profit;
    }

    void setTime(float time) {
        TOP::time = time;
    }

    void addTime(float time){
        TOP::time += time;
    }

    float distance_time(Vertex current, Vertex next){
        int x1 = current.getX();
        int y1 = current.getY();
        int x2 = next.getX();
        int y2 = next.getY();

        return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
    }

    void calculate_solution(vector<int>solution , int check){
        profit = 0;
        int length = solution.size();
        int visitedNodes = 0;
        for (int index = 0; index < length; index++){

//            cout<<time<<" <--Time"<<'\n';
//            cout<<profit<<"<--Profit"<<'\n';

            if(solution[index] < 0){
                setTime(0); // New path
                continue;
            }
            else{
                Vertex current = vertexVector[solution[index]]; //Current Vertex
                int canBeVisitedResult = canBeVisited(current);
                if(canBeVisitedResult == 0){
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
//                    cout<<"Can visit "<<current.getI()<<'\n';
                    float duration = current.getD();
                    addTime(duration);
                    addProfit(current.getProfit());

                }
                else if(canBeVisitedResult == -1){
//                    cout<<"Cant's visit node "<<current.getI()<<'\n';
                }

                else{
                    // Wait for a time
//                    cout<<"Can visit "<<current.getI()<<" But after waiting "<<canBeVisitedResult<<'\n';
                    addTime(canBeVisitedResult);
                    float duration = current.getD();
                    addTime(duration);
                    addProfit(current.getProfit());
                    if(solution[index] != 0){
                        visitedNodes++;
                    }
                }
                if(index + 1 != length){
                    Vertex next = vertexVector[solution[index + 1]];
                    float travelTime = distance_time(current, next);
                    addTime(travelTime);
                }

            }
            }

        if(check == 1){
            cout<<"We visit "<<visitedNodes<<" nodes with profit "<<profit<<'\n';
        }

    }

    int canBeVisited(Vertex v){
        if (time >= v.getOpeningTime() && time <= v.getClosingTime())
            return 0;
        else if (time < v.getOpeningTime())
            return v.getOpeningTime() - time;
        else{
            return -1;
        }
    }

    vector<int> random_solution_generator(){
        vector<int> numbers;
        unsigned num = chrono::system_clock::now().time_since_epoch().count();
        for(int i = 1; i < N+1; i++)
            numbers.push_back(i);

        shuffle(numbers.begin(), numbers.end(),default_random_engine(num));

        vector<int> solution;
        solution.push_back(0); // We start from 0
        solution.push_back(numbers[0]);
        float nTest = N;
        float vTest = V;
        int division = ceil(nTest/vTest);
        for (int i = 1; i<N; i++){
            if ((i % division)==0){
                solution.push_back(0); //Each path ends to 0
                solution.push_back(-1); // To separate path
                solution.push_back(0); // next path start
            }
            solution.push_back(numbers[i]);
        }
        solution.push_back(0); //Last part
        return solution;

    }

    vector<int> grasp_solution_generator() {
        vector<int> solution;
        solution.push_back(0); // We start from 0
        solution.push_back(0);
        float nTest = N;
        float vTest = V;
        int division = ceil(nTest/vTest);
        int pathN = -1;
        for (int i = 1; i<N; i++){
            if ((i % division)==0){
                solution.push_back(0); //Each path ends to 0
                solution.push_back(pathN); // To separate path
                solution.push_back(0); // next path start
                pathN --;
            }
            solution.push_back(0);
        }
        solution.push_back(0); //Last part
        return solution;
    }
};

class SA{
protected:
    int maxIterations;
    float t0;
    float alpha;
    float noImprove;
    int B;
    int N;
    int V;
    vector<int> solution;
    int fBest;

public:

    SA(vector<int> solution_ , float t0_ , float  alpha_ , float noImprove_, int B_, int N_, int V_){
        solution = solution_;
        alpha = alpha_;
        t0 = t0_;
        B = B_;
        N = N_;
        V = V_;
        fBest = 0;
        noImprove = noImprove_;
        maxIterations = (N_ + V_ -1) * B_;
    }

    vector<int> swapRandom(vector<int> mySolution){

        vector<int> newSolution = mySolution;
        int randomIndex1 = vertexGenerator(mySolution);
        int randomIndex2 = vertexGenerator(mySolution);
        int tmp = newSolution[randomIndex1];
        newSolution[randomIndex1] = newSolution[randomIndex2];
        newSolution[randomIndex2] = tmp;
        return newSolution;
    }

    bool checkVertex(int v1){
        if(v1 == 0 || v1 == -1){
            return false;
        }
        return true;
    }

    int vertexGenerator(vector<int> mySolution){

        int length = mySolution.size();
        int randomIndex = (rand() % length);
        int vertex = mySolution[randomIndex];
        while (!checkVertex(vertex)){
            randomIndex = (rand() % length);
            vertex = mySolution[randomIndex];
        }
        return randomIndex;

    }

    vector<int> insertionRandom(vector<int> mySolution){
        vector<int> newSolution = mySolution;
        int randomIndex1 = vertexGenerator(mySolution);
        int randomIndex2 = vertexGenerator(mySolution);
        newSolution.erase(newSolution.begin() + randomIndex1);
        if(randomIndex2 > randomIndex1){
            auto position = newSolution.begin() + randomIndex2 - 1;
            newSolution.insert(position , mySolution[randomIndex1]);
        }
        else{
            auto position = newSolution.begin() + randomIndex2;
            newSolution.insert(position , mySolution[randomIndex1]);
        }
        return newSolution;

    }

    int random_0_1(){
        return (rand()%100);
    }

    float calculate(vector<int> someSolution , int check){
        TOP top(0);
        top.calculate_solution(someSolution, check);
        return top.getProfit();
    }

    vector<int> swapVector(vector<int> vector, int index1, int index2){
        int tmp = vector[index1];
        vector[index1] = vector[index2];
        vector[index2] = tmp;
        return vector;
    }

    vector<int> insertionVector(vector<int> vector1, int index1 , int index2){
        vector<int> newVector = vector1;
        newVector.erase(newVector.begin() + index1);
        if(index2 > index1){
            auto position = newVector.begin() + index2 - 1;
            newVector.insert(position , vector1[index1]);
        }
        else{
            auto position = newVector.begin() + index2;
            newVector.insert(position , vector1[index1]);
        }
        return newVector;
    }

    bool localSearchSwap(){
        vector<int> tempSolution = solution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate(solution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] == 0 || tempSolution[i] == -1){ //We don't want to swap these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] == 0 || tempSolution[j] == -1){
                    continue;
                }
                tempSolution = swapVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = solution;
            }
        }

        if (flag == 1) {
            solution = swapVector(solution, bestI, bestJ);

            return true;
        }

        return false;

    }

    bool localSearchInsertion(){
        vector<int> tempSolution = solution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate(solution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] == 0 || tempSolution[i] == -1){ //We don't want to mess with these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] == 0 ||  tempSolution[j] == -1){
                    continue;
                }
                tempSolution = insertionVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = solution;
            }
        }

        if (flag == 1) {
            solution = insertionVector(solution, bestI, bestJ);
            return true;
        }

        return false;
    }

    void ssaAlgorithm(){
        // Initial Part!
        srand(time(0));
        vector<int> currentSolution = solution;
        float t = t0;
        fBest = calculate(solution , 1);
        int localNoImprove = 0;
        // Initial Part!
        vector<int> newSolution; // Y
        while(localNoImprove < noImprove){
            int iteration = 0;
            while(iteration < maxIterations){

                iteration++; // Inc iteration
                float currentProfit = calculate(currentSolution , 0); // OBJ(X)
                int p = random_0_1();
                if(p > 50){
                    newSolution = insertionRandom(currentSolution);
                }
                else{
                    newSolution = swapRandom(currentSolution);
                }

                float newProfit = calculate(newSolution , 0);
                float delta = newProfit - currentProfit;
                if(delta > 0){
                }
                else{
                    double  x = delta / t;
                    double prob = exp(x) * 100;
                    float r = random_0_1();
                    if (prob < r){
                        continue;
                    }
                }
                currentSolution.clear();
                currentSolution = newSolution;

                if(newProfit > fBest){
                    solution = currentSolution;
                    fBest = newProfit;
                }
                newSolution.clear();

            }
            t = t * alpha;
            // Local search
            bool isSwapImprove = localSearchSwap();
            bool isInsertionImprove = localSearchInsertion();
            if(!isSwapImprove && !isInsertionImprove){
                localNoImprove++;
            }
        }

    }

    float printFinalSolution(){
        cout<<"***************** Final Solution *****************"<<'\n';
        for(int k = 0; k<solution.size(); k++){
            cout<<solution[k]<<' ';
        }
        cout<<'\n';

        float finalProfit = calculate(solution, 1);
        return finalProfit;
    }

    void fsaAlgorithm(int seconds){
        // Initial Part!
        srand(time(0));
        vector<int> currentSolution = solution;
        float t = t0;
        fBest = calculate(solution , 1);
//        int localNoImprove = 0;
        // Initial Part!
        vector<int> newSolution; // Y
        auto Start = std::chrono::high_resolution_clock::now();
        while(1){
            int iteration = 0;
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> Elapsed = End - Start;
            if (Elapsed.count() >= seconds * 1000)
                break;
            while(iteration < maxIterations){

                iteration++; // Inc iteration
                float currentProfit = calculate(currentSolution , 0); // OBJ(X)
                int p = random_0_1();
                if(p > 50){
                    newSolution = insertionRandom(currentSolution);
                }
                else{
                    newSolution = swapRandom(currentSolution);
                }

                float newProfit = calculate(newSolution , 0);
                float delta = newProfit - currentProfit;
                if(delta > 0){
                }
                else{
                    double  x = delta / t;
                    double prob = exp(x) * 100;
                    float r = random_0_1();
                    if (prob < r){
                        continue;
                    }
                }
                currentSolution.clear();
                currentSolution = newSolution;

                if(newProfit > fBest){
                    solution = currentSolution;
                    fBest = newProfit;
                }
                newSolution.clear();

            }
            t = t * alpha;
            // Local search
            bool isSwapImprove = localSearchSwap();
            bool isInsertionImprove = localSearchInsertion();
        }

    }

};

class GRASP{
protected:
    int N;
    int V;
    int RCL_LEN;
    int maxIteration;
    vector<int> FirstSolution;
    vector<int> mainSolution;
    vector<int> BestSolution;
    int maxLocalI;
    int fBest;
public:
    GRASP(int n_, int v_, int rcl_len, int maxI, int maxLI, vector<int> FS){
        N = n_;
        V = v_;
        RCL_LEN = rcl_len;
        maxIteration =  maxI;
        fBest = 0;
        FirstSolution = FS;
        maxLocalI = maxLI;
    }

    struct sortByProfitStructure
    {
        inline bool operator() (const Vertex& v1, const Vertex& v2)
        {
            return (v1.getProfit() > v2.getProfit());
        }
    };

    struct sortByTravelTimeStructure{
        inline bool operator() (const Vertex& v1, const Vertex& v2){
            TOP top(0);
            Vertex initial = vertexVector[0];
            float v1TravelTime = top.distance_time(initial, v1);
            float v2TravelTime = top.distance_time(initial, v2);
            return (v1TravelTime < v2TravelTime);
        }
    };

    vector<int> sortByProfit(){
        vector<Vertex> tmpV = vertexVector;
        sort(tmpV.begin(), tmpV.end(), sortByProfitStructure());
        vector<int> sortedByProfit;
        for(int i = 0; i < tmpV.size() - 1; i++){
            sortedByProfit.push_back(tmpV[i].getI());
        }
        return sortedByProfit;
    }

    vector<int> sortByTravelTime(){
        vector<Vertex> tmpV = vertexVector;
        sort(tmpV.begin(), tmpV.end(), sortByTravelTimeStructure());
        TOP top(0);
        vector<int> sortedByTT;
        for(int i = 1; i< tmpV.size(); i++){
            sortedByTT.push_back(tmpV[i].getI());
        }

        return sortedByTT;
    }

   int randomZeroRCL(int size){
        int randomNumber = rand() % size;
       return randomNumber;
    }

    int randomZero100(){
        int randomNumber = rand() % 100;
        return randomNumber;
    }

    vector<int> addToRCL( vector<int> sortedV){
        vector<int> RCL;
        int index = 0;
        for(int & vertexId : sortedV){
            if(index < RCL_LEN){
                RCL.push_back(vertexId);
                index++;
            }
            else{
                break;
            }
        }
        return RCL;
    }

    float calculate(vector<int> someSolution , int check){
        TOP top(0);
        top.calculate_solution(someSolution, check);
        return top.getProfit();
    }

    void constructSolution(){
        vector<int> greedyVector;
        mainSolution = FirstSolution;
        int chooseGreedy = randomZero100();
        if (chooseGreedy > 50){
            greedyVector = sortByTravelTime();
        } else{
            greedyVector = sortByProfit();
        }
        for(int i = 0; i < N; i++){ //Construct Path
            vector<int> RCL = addToRCL(greedyVector);
            int index = randomZeroRCL(RCL.size());
            int selectedCandidate = RCL[index];
            mainSolution = choosingBestPath(mainSolution, selectedCandidate);
            RCL.erase(RCL.begin() + index); //Delete from RCL
            greedyVector.erase(greedyVector.begin() + index); //Delete From vector
        }
//        printVector("Constructed Solution", mainSolution);
//        calculate(mainSolution, 1);
    }

    void printVector(string name, vector<int> v){
        cout<<name<<'\n';
        for(int & element: v){
            cout<<element<<" ";
        }
        cout<<'\n';
    }

    vector<int> choosingBestPath(vector<int> baseSolution, int SC){
        float bestProfit = calculate(baseSolution,0);
        int bestI = 0;
        int ableToadI = 0;
        int flag = 0;
        for(int i = 0 ; i< V; i++){
            vector<int> tmp = baseSolution;
            int startIndex = 1 + ((i) * 3) + (i * ceill(N/V));
            int lastIndex = startIndex + ceill(N/V);
            for (int j = startIndex ; j < lastIndex ; j++){
                if(tmp[j] == 0){
                    tmp[j] = SC;
                    ableToadI = i;
                    break;
                }
                continue;
            }
            float currentProfit = 0;
            currentProfit = calculate(tmp, 0);
            if(currentProfit >  bestProfit){
                bestProfit = currentProfit;
                bestI = i;
                flag = 1;
            }
            else{
                if(flag == 0)
                    bestI = ableToadI;
            }
        }

        // Now add to the best position
        int startIndexFinal = 1 + ((bestI) * 3) + (bestI * ceill(N/V));
        int lastIndexFinal = startIndexFinal + ceill(N/V);
        for (int j = startIndexFinal ; j < lastIndexFinal; j++){
            if(baseSolution[j] == 0){
                baseSolution[j] = SC;
                break;
            }
        }
        return baseSolution;
    }

    vector<int> swapVector(vector<int> vector, int index1, int index2){
        int tmp = vector[index1];
        vector[index1] = vector[index2];
        vector[index2] = tmp;
        return vector;
    }

    vector<int> insertionVector(vector<int> vector1, int index1 , int index2){
        vector<int> newVector = vector1;
        newVector.erase(newVector.begin() + index1);
        if(index2 > index1){
            auto position = newVector.begin() + index2 - 1;
            newVector.insert(position , vector1[index1]);
        }
        else{
            auto position = newVector.begin() + index2;
            newVector.insert(position , vector1[index1]);
        }
        return newVector;
    }

    bool localSearchSwap(){
        vector<int> tempSolution = mainSolution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate(mainSolution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 1){ //We don't want to swap these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 1){
                    continue;
                }
                tempSolution = swapVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = mainSolution;
            }
        }

        if (flag == 1) {
            mainSolution = swapVector(mainSolution, bestI, bestJ);

            return true;
        }

        return false;

    }

    bool localSearchInsertion(){
        vector<int> tempSolution = mainSolution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate(mainSolution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 1){ //We don't want to mess with these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 1){
                    continue;
                }
                tempSolution = insertionVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = mainSolution;
            }
        }

        if (flag == 1) {
            mainSolution = insertionVector(mainSolution, bestI, bestJ);
            return true;
        }

        return false;
    }

    void graspAlgorithm(){
        srand(time(0));
        BestSolution = FirstSolution;
        float bestProfit = calculate(BestSolution, 0);
        for(int j = 0; j< maxIteration; j++){
            constructSolution();
            for(int i = 0; i < maxLocalI; i++){
                localSearchSwap();
                localSearchInsertion();
            }
            float currentProfit = calculate(mainSolution, 0);
            if(currentProfit > bestProfit){
                BestSolution = mainSolution;
                bestProfit = currentProfit;
            }
        }
        printVector("Final Solution" , BestSolution);
        calculate(BestSolution, 1);
    }

    void graspAlgorithmTime(int seconds){
        srand(time(0));
        BestSolution = FirstSolution;
        float bestProfit = calculate(BestSolution, 0);
        auto Start = std::chrono::high_resolution_clock::now();
        while(1){
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> Elapsed = End - Start;
            if (Elapsed.count() >= seconds * 1000)
                break;
            constructSolution();
            for(int i = 0; i < maxLocalI; i++){
                localSearchSwap();
                localSearchInsertion();
            }
            float currentProfit = calculate(mainSolution, 0);
            if(currentProfit > bestProfit){
                BestSolution = mainSolution;
                bestProfit = currentProfit;
            }
        }
        printVector("Final Solution" , BestSolution);
        calculate(BestSolution, 1);
    }

};

class TabuItem{
private:
    int I;
    int J;
    float changeProfit;
    float actualProfit;
public:
    TabuItem(int i, int j, float changeProfit, float actualProfit) : I(i), J(j), changeProfit(changeProfit),
                                                                     actualProfit(actualProfit) {}

    int getI() const {
        return I;
    }

    int getJ() const {
        return J;
    }

    float getChangeProfit() const {
        return changeProfit;
    }

    float getActualProfit() const {
        return actualProfit;
    }

};

class TabuSearch {
protected:
    int Tabu_Len;
    vector<int> FirstSolution;
    vector<int> MainSolution;
    vector<int> FinalSolution;
    int maxIteration;
    int fBest;
    vector<TabuItem> tabuList;
    vector<TabuItem> neighborhood;
public:
    TabuSearch (int tabu_len, vector<int> fs, int maxI){
        Tabu_Len = tabu_len;
        FirstSolution = fs;
        maxIteration = maxI;
        fBest = 0;
    }

    float calculate_profit(){
        int profitSum = 0;
        for(int i = 0 ; i< vertexVector.size(); i++){
            profitSum = profitSum + vertexVector[i].getProfit();
        }
        return profitSum;
    }

    bool checkTabuList(TabuItem tabu){
        for(int i = 0; i< tabuList.size(); i++){
            if(tabuList[i].getI() == tabu.getI() && tabuList[i].getJ() == tabu.getJ()){
                return false; // Item is in tabu list
            }
            else if(tabuList[i].getI() == tabu.getJ() && tabuList[i].getJ() == tabu.getI()){
                return false;
            }
        }
        return true;
    }

    float calculate(vector<int> someSolution , int check){
        TOP top(0);
        top.calculate_solution(someSolution, check);
        return top.getProfit();
    }

    struct sortByChangeStructure
    {
        inline bool operator() (const TabuItem& tb1, const TabuItem& tb2)
        {
            return (tb1.getChangeProfit() > tb2.getChangeProfit());
        }
    };

    vector<int> swapVector(vector<int> vector, int index1, int index2){
        int tmp = vector[index1];
        vector[index1] = vector[index2];
        vector[index2] = tmp;
        return vector;
    }

    vector<int> insertionVector(vector<int> vector1, int index1 , int index2){
        vector<int> newVector = vector1;
        newVector.erase(newVector.begin() + index1);
        if(index2 > index1){
            auto position = newVector.begin() + index2 - 1;
            newVector.insert(position , vector1[index1]);
        }
        else{
            auto position = newVector.begin() + index2;
            newVector.insert(position , vector1[index1]);
        }
        return newVector;
    }

    int localSearchSwap(){
        vector<int> tempSolution = MainSolution;
        int solutionLength = tempSolution.size();
        float currentProfit = calculate(MainSolution , 0);
        neighborhood.clear(); // Clear all data from neighborhood;
        float totalProfit = calculate_profit();

        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 1){ //We don't want to swap these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 1){
                    continue;
                }
                tempSolution = swapVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                float profitChange = localProfit - currentProfit;
                createMoves(i,j,profitChange, localProfit); // Add this move to the list
                tempSolution.clear();
                tempSolution = MainSolution;
            }
        }

        sort(neighborhood.begin(),neighborhood.end(),sortByChangeStructure()); //Sort them by profit
        for(int index = 0; index < neighborhood.size(); index++){
            if(checkTabuList(neighborhood[index])) { //Item is not in tabu list
                // We can use that
                // Add it to tabu list first
                addToTabuList(neighborhood[index]);
                int I = neighborhood[index].getI();
                int J = neighborhood[index].getJ();
                MainSolution = swapVector(MainSolution,I,J); // We change solution to them!
                return 0;
                break;
            }
            else{
                if(neighborhood[index].getActualProfit() > 0.98 * totalProfit ){
                    addToTabuList(neighborhood[index]);
                    int I = neighborhood[index].getI();
                    int J = neighborhood[index].getJ();
                    MainSolution = swapVector(MainSolution,I,J); // We change solution to them!
                    return 1;
                }
                continue; //Item is in tabu list
            }
        }
    }

    bool localSearchInsertion(){
        vector<int> tempSolution = MainSolution;
        int solutionLength = tempSolution.size();
        float bestProfit = calculate(MainSolution , 0);
        int bestI;
        int bestJ;
        int flag = 0;
        for(int i = 0; i < solutionLength; i++){
            if(tempSolution[i] < 1){ //We don't want to mess with these
                continue;
            }
            for(int j = i; j < solutionLength; j++){
                if(tempSolution[j] < 1){
                    continue;
                }
                tempSolution = insertionVector(tempSolution, i, j);
                float localProfit = calculate(tempSolution , 0);
                if(localProfit > bestProfit){
                    flag = 1;
                    bestProfit = localProfit;
                    bestI = i;
                    bestJ = j;
                }
                tempSolution.clear();
                tempSolution = MainSolution;
            }
        }

        if (flag == 1) {
            MainSolution = insertionVector(MainSolution, bestI, bestJ);
            return true;
        }

        return false;
    }

    void createMoves(int i, int j, float profitChange, float actualProfit){
        TabuItem tabuItem(i,j,profitChange, actualProfit);
        neighborhood.push_back(tabuItem);
    }

    void addToTabuList(TabuItem tabu){
        if(tabuList.size() == Tabu_Len){
            tabuList.erase(tabuList.begin()); //We remove first element
        }
        tabuList.push_back(tabu);
    }

    void printVector(string name, vector<int> v){
        cout<<name<<'\n';
        for(int & element: v){
            cout<<element<<" ";
        }
        cout<<'\n';
    }

    void tabuSearchAlgorithm(){
        float firstProfit = calculate(FirstSolution,1);
        printVector("First Solution", FirstSolution);
        MainSolution = FirstSolution;
        for(int i = 0; i < maxIteration; i++){
            int response = localSearchSwap(); //Main tabu search operator
            if(response == 1){
                //We have something greater than 97 percent
                break;
            }
        }

        bool upgradeResult = localSearchInsertion();
        while(upgradeResult){
            upgradeResult = localSearchInsertion();
        }
        printVector("Final Solution", MainSolution);
        calculate(MainSolution, 1);
    }

    void tabuSearchAlgorithmTime(int seconds){
        float firstProfit = calculate(FirstSolution,1);
        printVector("First Solution", FirstSolution);
        MainSolution = FirstSolution;
        FirstSolution = FirstSolution;
        fBest = firstProfit;
        auto Start = std::chrono::high_resolution_clock::now();
        while (1) {
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> Elapsed = End - Start;
            if (Elapsed.count() >= seconds * 1000)
                break;
            int response = localSearchSwap(); //Main tabu search operator
            if (response == 1) {
                //We have something greater than 97 percent
                break;
            }
            bool upgradeResult = localSearchInsertion();
            while (upgradeResult) {
                upgradeResult = localSearchInsertion();
            }

            float currentProfit = calculate(MainSolution,0);
            if(currentProfit > fBest){
                FinalSolution = MainSolution;
                fBest = currentProfit;
            }

        }
        printVector("Final Solution", FinalSolution);
        calculate(FinalSolution, 1);
    }


};

void add_to_vertex_vector(Vertex v){
    vertexVector.push_back(v);
}

float calculate_profit(){
    int profitSum = 0;
    for(int i = 0 ; i< vertexVector.size(); i++){
        profitSum = profitSum + vertexVector[i].getProfit();
    }
    return profitSum;
}

#endif //HW_1_INTERFACE_H
