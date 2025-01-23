#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace seal;
using namespace std;


int main() {


    int ring_dim = 32768;
    int p = 65537;

    vector<int> x_vec(32768, 2);
    vector<int> coeff_vec(32768, 18);

    vector<vector<int>> x_vec_power(512, vector<int>(32768));

    x_vec_power[0] = x_vec;

    // vector<int> x_vec(65536, 0);
    // vector<int> coeff(65536, 18);

    chrono::high_resolution_clock::time_point time_start, time_end;
    time_start = chrono::high_resolution_clock::now();

    // for (int i = 1; i < (int) x_vec.size(); i++) {
    //     x_vec[i] = (x_vec[0] * x_vec[i-1]) % p;
    // }

    // for (int i = 0; i < (int) x_vec.size(); i++) {
    //     x_vec[i] = (x_vec[i] * coeff[i]) % p;
    // }

    // for (int i = 1; i < (int) x_vec.size(); i++) {
    //     x_vec[0] = (x_vec[i] + x_vec[0]) % p;
    // }


    for (int i = 1; i < 512; i++) {
        for (int j = 0; j < ring_dim; j++) {
            x_vec_power[i][j] = (x_vec_power[i-1][j] * x_vec_power[i-1][j]) % p;
        }
    }

    for(int i = 0; i < 256; i++) {
        for(int j = 0; j < 256; j++) {
            for (int c = 0; c < ring_dim; c++) {
                x_vec_power[j+256][c] += (coeff_vec[c] * x_vec_power[j][c]) % p; // one multi + one add
            }
        }

        for (int c = 0; c < ring_dim; c++) {
            x_vec_power[i][c] += (coeff_vec[c] * x_vec_power[i][c]) % p; // one multi + one add
        }
    }

    time_end = chrono::high_resolution_clock::now();

    cout << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;




}