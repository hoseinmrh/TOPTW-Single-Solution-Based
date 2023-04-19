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

    int getN() const {
        return N;
    }

    int getV() const {
        return V;
    }

    float getTime() const {
        return time;
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

            if(solution[index] == -1){
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
};

class SA{
protected:
    int maxIterations;
    float t0;
    float alpha;
    float noImprove;
    float maxT;
    int B;
    int N;
    int V;
    vector<int> solution;
    int fBest;

public:
    SA(vector<int> solution_ ){
        solution = solution_;
    }

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

    const vector<int> &getSolution() const {
        return solution;
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
