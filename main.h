//
// Created by peterchaozhu on 16/2/21.
//

#ifndef OPTION_MAIN_H
#define OPTION_MAIN_H
#include <vector>
#include<iostream>
#include <math.h>

using namespace std;
class Normality{
    /* normality stats class for data normality test */
public:
    // data
    float mean;
    float sum;
    long count;
    float var;
    float skew;
    float kurt;
    float JB_factor;
    // methods
    float get_sum(vector<float>& return_v);
    float get_mean();
    void update(vector<float>& return_v);
    float get_JB_factor();
    // constructors
    Normality();
    Normality(vector<float>& return_v);
};
float Normality::get_sum(vector<float>& vec) {
    /* derive the sum of normality */
    while (count < vec.size()) {
        sum = sum + vec[count];
        count ++;
    }
    return sum;
}
float Normality::get_mean() {
    /* derive the mean from sum */
    mean = sum / (float)count;
    return mean;
}
void Normality::update(vector<float>& vec)
{
    /* update stats info inside the class */
    for (int i = 0; i < count; i++) {
        float temp_std = vec[i] - mean;
        var = pow(temp_std, 2.0) + var;
        skew = pow(temp_std, 3.0) + skew;
        kurt = pow(temp_std, 4.0) + kurt;
    }
    var = var / count;
    skew = skew / count;
    kurt = kurt / count;
    skew = skew/ (pow (var, 1.5));
    kurt = kurt/ (pow (var, 2.0));
}
float Normality::get_JB_factor() {
    /* calculate JB factor */
    JB_factor = kurt - 3.0;
    JB_factor = pow(JB_factor, 2.0);
    JB_factor = JB_factor / 4.0;
    JB_factor = JB_factor + pow(skew, 2.0);
    JB_factor = JB_factor * count;
    JB_factor = JB_factor / 6.0;
    return JB_factor;
}
Normality::Normality(){
    /* initialization constructor */
    mean = 0.0;
    sum = 0.0;
    count = 0.0;
    var = 0.0;
    skew = 0.0;
    kurt = 0.0;
    JB_factor = 0.0;
}
Normality::Normality(vector<float>& vec){
    /* complicated constructor */
    sum = get_sum(vec);
    mean = get_mean();
    update(vec);
    JB_factor = get_JB_factor();
}

#endif //OPTION_MAIN_H
