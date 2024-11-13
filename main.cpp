#include"compcore.hpp"
#include <ctime>
#include <random>
#include <chrono>
#include<iostream>

using namespace std;

typedef vector<float> VectorDouble;
typedef vector<VectorDouble> MatrixDouble;
typedef vector<int> VectorInt;
typedef vector<VectorInt> MatrixInt;

class _LSA_Data
{
    public:
    int max_shift;
    VectorDouble X;
    VectorDouble Y;
    _LSA_Data(){ VectorDouble X; VectorDouble Y; max_shift=numeric_limits<int>::infinity(); };
    void _assign(int, VectorDouble, VectorDouble);
};

void _LSA_Data::_assign(int shift, VectorDouble x, VectorDouble y)
{
    max_shift = shift;
    X.assign(x.begin(),x.end());
    Y.assign(y.begin(),y.end());
}

class LSA_Result
{
    public:
      float score;
      MatrixInt trace;
};

LSA_Result DP_lsa( const _LSA_Data& data, bool keep_trace )
{
    LSA_Result lsa_result;
    int max_p[2] = { 0 }; 
    int porn = 0;
    float max_s = -numeric_limits<float>::infinity();

    MatrixDouble psm = vector<vector<float> >(data.X.size() + 1, vector<float>(data.Y.size() + 1));
    MatrixDouble nsm = vector<vector<float> >(data.X.size() + 1, vector<float>(data.Y.size() + 1));

    for (int i = 1; i <= data.X.size(); i++)
        for (int j = max(1, i - data.max_shift); j <= min((int)data.Y.size(), (int)i + data.max_shift); j++)
        {
            float s1 = data.X[i - 1] * data.Y[j - 1];
            psm[i][j] = max(float(0.0), psm[i - 1][j - 1] + s1);
            nsm[i][j] = max(float(0.0), nsm[i - 1][j - 1] - s1);

            if (psm[i][j] >= max_s)
            {
                max_p[0] = i; max_p[1] = j; max_s = psm[i][j]; porn = 1;
            }
            if (nsm[i][j] >= max_s)
            {
                max_p[0] = i; max_p[1] = j; max_s = nsm[i][j]; porn = -1;
            }
        }

    int length = 0; 
    vector<int> step; 
    step.push_back(max_p[0]);
    step.push_back(max_p[1]);

    if (porn == -1)
    {
        lsa_result.score = -1 * nsm[max_p[0]][max_p[1]] / data.X.size();
        while (nsm[max_p[0] - length][max_p[1] - length] != 0. && keep_trace == true)
        {
            length++; lsa_result.trace.push_back(step); step.clear();
            step.push_back(max_p[0] - length); step.push_back(max_p[1] - length);
        }
    }
    else
    {
        lsa_result.score = psm[max_p[0]][max_p[1]] / data.X.size();
        while (psm[max_p[0] - length][max_p[1] - length] != 0. && keep_trace == true)
        {
            length++; lsa_result.trace.push_back(step); step.clear();
            step.push_back(max_p[0] - length); step.push_back(max_p[1] - length);
        }
    }
    return lsa_result;
}


// #define COL 10 // 20, 50，100， 200， 500， 1000 (gpu 7次,cpu7次)
#define shift 10

int main()
{
    vector<double>new_algorithm;
    vector<double>old_algorithm;
    vector<int>N;
    vector<int>col;

    col.push_back(10);
    col.push_back(20);
    col.push_back(50);
    col.push_back(100);
    col.push_back(200);
    col.push_back(500);
    col.push_back(1000);


    N.push_back(2000);

    for(const auto& COL : col){


        for (const auto& element : N) 
        {
            int num = element;

            vector<float>a_num(num * COL);
            default_random_engine e;
            uniform_real_distribution<float> u(-1,1);
            e.seed(time(0));
            for (int i = 0; i < num * COL; ++i)
                a_num[i] = u(e);

            auto start_time_0 = chrono::high_resolution_clock::now();
            vector<float>X(a_num);
            vector<float>Y=(X);



            LSA lsa;
            lsa.assign(num, num, shift, COL, X, Y);
            lsa.dp_lsa();
            lsa.lsa_clean();

            auto end_time_0 = chrono::high_resolution_clock::now();
            auto duration_0 = chrono::duration_cast<chrono::microseconds>(end_time_0 - start_time_0);
            

            auto start_time = chrono::high_resolution_clock::now();
            _LSA_Data Data;
            LSA_Result res;
            vector<float>res_score;
            vector<int>x_0;
            vector<int>x_1;
            vector<int>y_0;
            vector<int>y_1;
            for (int i = 0; i < num; ++i)
            {
                vector<float>X;
                X.assign(a_num.begin()+ i * COL, a_num.begin() + (i + 1) * COL);
                for (int j = 0; j < num; ++j)
                {
                    vector<float>Y;


                    Y.assign(a_num.begin()+ j * COL, a_num.begin() + (j + 1) * COL);
                    Data._assign(shift, X, Y);

                    res = DP_lsa(Data, 1);

                    res_score.push_back(res.score);
                    x_0.push_back(res.trace.back()[0]);
                    x_1.push_back(res.trace.front()[0]);
                    y_0.push_back(res.trace.back()[1]);
                    y_1.push_back(res.trace.front()[1]);
                }
            }

            auto end_time = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

            // cout << "Time taken by new_design algorithm: " << duration_0.count() << " microseconds" << endl;
            // cout << "Time taken by old_design algorithm: " << duration.count() << " microseconds" << endl;

            new_algorithm.push_back(duration_0.count());
            old_algorithm.push_back(duration.count());
        }
    
    }
    
    cout << "                               num: ";
    for (const auto& element : col)
        cout << element << "\t";
        cout << "\n";

    cout << "Time taken by new_design algorithm: ";
    for (const auto& element : new_algorithm)
        cout << element << "\t";
        cout << "\n";

    cout << "Time taken by old_design algorithm: ";
    for (const auto& element : old_algorithm)
        cout << element << "\t";
        cout << "\n";

}

