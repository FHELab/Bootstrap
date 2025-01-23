#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>
#include <algorithm>

using namespace seal;
using namespace std;


// used to bootstrap for BFV ciphertexts, encrypting [0, t-1, r], where t is the prime, and r is the interval (≥ error bound = 128 in our case)
int main() {

    ////////////////////////////////////////////// PREPARE (R)LWE PARAMS ///////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 32768;
    BootstrapParam bootstrap_param = BootstrapParam(65537, 128, 512, 256, 256);
    int p = bootstrap_param.ciphertextSpacePrime;
    vector<uint64_t> rangeCheckIndices = rangeCheckIndices_bfv;

    int iter_cnt = 25;
    int group_cnt = 1;
    int skk_cnt = 1;

    vector<int> avg_success(6, 0);
    vector<int> max_success(6, 0);
    vector<int> min_success(6, n);
    vector<int> max_same_wrong(6, 0);
    vector<int> max_all_wrong(6, 0);
    vector<int> avg_same_wrong(6, 0);
    vector<int> avg_all_wrong(6, 0);

    // vector<int> threshold_list = {60, 72, 80, 96, 108, 120};

    int threshold_1 = 60;
    int threshold_2 = 72;
    int threshold_3 = 84;
    int threshold_4 = 96;
    int threshold_55 = 108;
    int threshold_6 = 120;


    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);

    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 60, 55 });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    ofstream datafile;
    datafile.open ("../ind.txt");

    // vector<float> threshold_5={-0.0002855605, -0.0002552205}; // for randomized rounding
    // vector<float> threshold_10={-0.0002822985, -0.0002484185}; // for randomized rounding
    // vector<float> threshold_15={-0.0002822985, -0.0002552205}; // for randomized rounding
    // vector<float> threshold_5 = {0.00426, -0.00626}; // ring_dim = 1024, threshold = 2
    // vector<float> threshold_10 = {0.00477, -0.00717}; // ring_dim = 1024, threshold = 3
    // vector<float> threshold_15 = {0.00582, -0.00822}; // ring_dim = 1024, threshold = 5
    vector<float> threshold_5 = {0.00091511, -0.00083722}; // ring_dim = 32768, threshold = 6
    vector<float> threshold_10 = {0.0011071714, -0.0010272686}; // ring_dim = 32768, threshold = 18
    vector<float> threshold_15 = {0.00124111075, -0.00116242425}; // ring_dim = 32768, threshold = 30
    vector<float> threshold_20 = {0.00152535295, -0.00144587205}; // ring_dim = 32768, threshold = 30
    vector<float> threshold_25 = {0.00175004445, -0.00166921055}; // ring_dim = 32768, threshold = 30
    vector<float> threshold_30 = {0.0019788607, -0.0019034243}; // ring_dim = 32768, threshold = 30

    // vector<float> threshold_5 = {0.0010414803, -0.0008484447}; // ring_dim = 16384, threshold = 6
    // vector<float> threshold_10 = {0.0012232206, -0.0010290444}; // ring_dim = 16384, threshold = 18
    // vector<float> threshold_15 = {0.001263616, -0.001068044}; // ring_dim = 16384, threshold = 30
    // vector<float> threshold_5 = {0.00582, -0.00548};
    // vector<float> threshold_10 = {0.00587, -0.00603};
    // vector<float> threshold_15 = {0.00765, -0.00708};

    chrono::high_resolution_clock::time_point time_start, time_end;


    for (int kk = 0; kk < skk_cnt; kk++) {

        prng_seed_type seed_lwe;
        for (auto &i : seed_lwe) {
            i = random_uint64();
        }
        auto rng_lwe = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed_lwe));
        bfv_params.set_random_generator(rng_lwe);

        SEALContext seal_context_lwe(bfv_params, true, sec_level_type::none);

        KeyGenerator lwe_key_keygen(seal_context_lwe, n);
        SecretKey lwe_key = lwe_key_keygen.secret_key();

        vector<vector<int>> threshold_5_sk(group_cnt, vector<int>(n));
        vector<vector<int>> threshold_10_sk(group_cnt, vector<int>(n));
        vector<vector<int>> threshold_15_sk(group_cnt, vector<int>(n));
        vector<vector<int>> threshold_20_sk(group_cnt, vector<int>(n));
        vector<vector<int>> threshold_25_sk(group_cnt, vector<int>(n));
        vector<vector<int>> threshold_30_sk(group_cnt, vector<int>(n));

        vector<int> sk_5_cnt(group_cnt, 0);
        vector<int> sk_10_cnt(group_cnt, 0);
        vector<int> sk_15_cnt(group_cnt, 0);
        vector<int> sk_20_cnt(group_cnt, 0);
        vector<int> sk_25_cnt(group_cnt, 0);
        vector<int> sk_30_cnt(group_cnt, 0);

        inverse_ntt_negacyclic_harvey(lwe_key.data().data(), seal_context_lwe.key_context_data()->small_ntt_tables()[0]);
        auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
        auto lwe_sk = regevGenerateSecretKey(lwe_params);
        for (int i = 0; i < n; i++) {
            lwe_sk[i] = (uint64_t) lwe_key.data()[i] > (uint64_t) p ? p-1 : lwe_key.data()[i];
            if ((uint64_t) lwe_key.data()[i] > (uint64_t) p) big_prime_global = (uint64_t) lwe_key.data()[i]+1;
            // cout << lwe_key.data()[i] << " ";
        }
        // cout << "------------------------------\n" << lwe_sk << "\n------------------------------\n";

        seal::util::RNSIter new_key_rns(lwe_key.data().data(), ring_dim);
        ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context_lwe.key_context_data()->small_ntt_tables());

        for (int jj = 0; jj < group_cnt; jj++) {


        
            vector<double> threshold_5_a_avg(n, 0);
            vector<double> threshold_10_a_avg(n, 0);
            vector<double> threshold_15_a_avg(n, 0);
            vector<double> threshold_20_a_avg(n, 0);
            vector<double> threshold_25_a_avg(n, 0);
            vector<double> threshold_30_a_avg(n, 0);


            vector<double> threshold_5_a_tot(n, 0);
            vector<double> threshold_10_a_tot(n, 0);
            vector<double> threshold_15_a_tot(n, 0);
            vector<double> threshold_20_a_tot(n, 0);
            vector<double> threshold_25_a_tot(n, 0);
            vector<double> threshold_30_a_tot(n, 0);

            vector<int> a_cnt(6, 0); // for threshold 5, 10, 15

            diff_global.resize(ring_dim);
            for (int ii = 0; ii < iter_cnt; ii++) {
                // cout << ii << endl;

                /////////////////////////////////////////////////// BOOTSTRAP //////////////////////////////////////////////////////

                prng_seed_type seed;
                for (auto &i : seed) {
                    i = random_uint64();
                }
                auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
                bfv_params.set_random_generator(rng);

                SEALContext seal_context(bfv_params, true, sec_level_type::none);
                // cout << "primitive root: " << seal_context.first_context_data()->plain_ntt_tables()->get_root() << endl;

                KeyGenerator keygen(seal_context);
                SecretKey bfv_secret_key = keygen.secret_key();
                seal::util::ConstPolyIter secret_key_bfv(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

                PublicKey bfv_public_key;
                keygen.create_public_key(bfv_public_key);

                Encryptor encryptor(seal_context, bfv_public_key);
                BatchEncoder batch_encoder(seal_context);
                Decryptor decryptor(seal_context, bfv_secret_key);



                // generate a key switching key based on key_before and secret_key
                KSwitchKeys ksk_to_lwe, ksk_to_bfv;
                
                seal::util::ConstPolyIter secret_key_lwe(lwe_key.data().data(), ring_dim, coeff_modulus.size());

                lwe_key_keygen.generate_kswitch_keys(secret_key_bfv, 1, static_cast<KSwitchKeys &>(ksk_to_lwe), false); // used to switch from secret_key_bfv to secret_key_lwe
                keygen.generate_kswitch_keys(secret_key_lwe, 1, static_cast<KSwitchKeys &>(ksk_to_bfv), false); // used to switch from secret_key_lwe to secret_key_bfv
                ksk_to_lwe.parms_id() = seal_context.key_parms_id();
                ksk_to_bfv.parms_id() = seal_context.key_parms_id();

                Decryptor decryptor_lwe(seal_context, lwe_key);

                Evaluator evaluator(seal_context);   


                MemoryPoolHandle my_pool = MemoryPoolHandle::New();

                // cout << "After param gen.\n";

                ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

                // one switching key for one lwe_sk
                // Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);
                Ciphertext sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);

                Ciphertext bfv_input;
                Plaintext pl;
                pl.resize(poly_modulus_degree_glb);
                pl.parms_id() = parms_id_zero;
                // vector<uint64_t> input_v(poly_modulus_degree_glb);
                for (int i = 0; i < (int) poly_modulus_degree_glb; i++) {
                    pl.data()[i] = 54;
                }
                // batch_xencoder.encode(input_v, pl);
                encryptor.encrypt(pl, bfv_input);
                // for (int i = 0; i < 10; i++) {
                //     cout << bfv_input.data(0)[i] << " ";
                // }
                // cout << endl;

                while(seal_context.last_parms_id() != bfv_input.parms_id()){
                    evaluator.mod_switch_to_next_inplace(bfv_input);
                }

                // for (int i = 0; i < (int) ring_dim; i++) {
                //     cout << bfv_input.data(0)[i] << " ";
                // }
                // cout << endl;

                Ciphertext copy_coeff = bfv_input;
                auto ct_in_iter = util::iter(copy_coeff);
                ct_in_iter += bfv_input.size() - 1;
                seal::util::set_zero_poly(ring_dim, 1, bfv_input.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level
                evaluator.switch_key_inplace(bfv_input, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk_to_lwe), 0, my_pool);


                // decryptor_lwe.decrypt(bfv_input, pl);
                // for (int i = 0; i < (int) ring_dim; i++) {
                // cout << pl.data()[i] << " ";
                // }
                // cout << endl << endl;
                vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(bfv_input, ring_dim, n);


                vector<int> ttmsg;
                time_start = chrono::high_resolution_clock::now();
                regevDec(ttmsg, lwe_ct_results, lwe_sk, lwe_params);
                time_end = chrono::high_resolution_clock::now();
                // cout << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;


                // cout << "heihei\n" << ttmsg << endl;

                // fix lwe sk, and do:
                // if some value > 59, > 64, > 69, record corresponding a vec, need diff in double, repeat 5 or 10 times, take average, 3*1024 avg, fill with null
                // 

                for (int mi = 0; mi < (int) ttmsg.size(); mi++) {
                    if (ttmsg[mi] > 10000) continue;

                    if (ttmsg[mi] > threshold_1) {
                        a_cnt[0]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_5_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                    if (ttmsg[mi] > threshold_2) {
                        a_cnt[1]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_10_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                    if (ttmsg[mi] > threshold_3) {
                        a_cnt[2]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_15_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                    if (ttmsg[mi] > threshold_4) {
                        a_cnt[3]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_20_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                    if (ttmsg[mi] > threshold_55) {
                        a_cnt[4]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_25_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                    if (ttmsg[mi] > threshold_6) {
                        a_cnt[5]+=1;
                        for (int ai = 0; ai < n; ai++) {
                            threshold_30_a_tot[ai] += diff_global[mi][ai];
                        }
                    }
                }
            }

            cout << endl << a_cnt << endl;

            vector<float> max_5(3, -5);
            vector<float> max_10(3, -5);
            vector<float> max_15(3, -5);
            vector<float> min_5(3, 5);
            vector<float> min_10(3, 5);
            vector<float> min_15(3, 5);
            vector<float> tot_5(3, 0);
            vector<float> tot_10(3, 0);
            vector<float> tot_15(3, 0);
            vector<float> max_20(3, -5);
            vector<float> max_25(3, -5);
            vector<float> max_30(3, -5);
            vector<float> min_20(3, 5);
            vector<float> min_25(3, 5);
            vector<float> min_30(3, 5);
            vector<float> tot_20(3, 0);
            vector<float> tot_25(3, 0);
            vector<float> tot_30(3, 0);
            vector<int> a_tot(3,0);

            for (int i = 0; i < n; i++) {
                threshold_5_a_avg[i] = threshold_5_a_tot[i] / a_cnt[0];
                threshold_10_a_avg[i] = threshold_10_a_tot[i] / a_cnt[1];
                threshold_15_a_avg[i] = threshold_15_a_tot[i] / a_cnt[2];
                threshold_20_a_avg[i] = threshold_20_a_tot[i] / a_cnt[3];
                threshold_25_a_avg[i] = threshold_25_a_tot[i] / a_cnt[4];
                threshold_30_a_avg[i] = threshold_30_a_tot[i] / a_cnt[5];
            
                threshold_5_sk[jj][i] = threshold_5_a_avg[i] > threshold_5[0] ? 65536 : threshold_5_a_avg[i] < threshold_5[1] ? 1 : 0;
                max_5[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_5_a_avg[i] > max_5[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_5_a_avg[i] : max_5[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_5[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_5_a_avg[i] < min_5[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_5_a_avg[i] : min_5[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_5[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_5_a_avg[i];
                a_tot[(lwe_sk[i].ConvertToInt()+1) % 65537] += 1;

                // if (threshold_5_sk[jj][i] != lwe_sk[i]) cout << "WRONG: " << threshold_5_a_avg[i] << "     " << threshold_5_sk[jj][i] << ", " << lwe_sk[i] << endl;
                
                threshold_10_sk[jj][i] = threshold_10_a_avg[i] > threshold_10[0] ? 65536 : threshold_10_a_avg[i] < threshold_10[1] ? 1 : 0;
                max_10[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_10_a_avg[i] > max_10[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_10_a_avg[i] : max_10[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_10[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_10_a_avg[i] < min_10[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_10_a_avg[i] : min_10[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_10[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_10_a_avg[i];
                
                threshold_15_sk[jj][i] = threshold_15_a_avg[i] > threshold_15[0] ? 65536 : threshold_15_a_avg[i] < threshold_15[1] ? 1 : 0;
                max_15[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_15_a_avg[i] > max_15[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_15_a_avg[i] : max_15[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_15[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_15_a_avg[i] < min_15[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_15_a_avg[i] : min_15[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_15[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_15_a_avg[i];

                threshold_20_sk[jj][i] = threshold_20_a_avg[i] > threshold_20[0] ? 65536 : threshold_20_a_avg[i] < threshold_20[1] ? 1 : 0;
                max_20[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_20_a_avg[i] > max_20[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_20_a_avg[i] : max_20[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_20[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_20_a_avg[i] < min_20[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_20_a_avg[i] : min_20[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_20[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_20_a_avg[i];
                
                threshold_25_sk[jj][i] = threshold_25_a_avg[i] > threshold_25[0] ? 65536 : threshold_25_a_avg[i] < threshold_25[1] ? 1 : 0;
                max_25[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_25_a_avg[i] > max_25[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_25_a_avg[i] : max_25[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_25[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_25_a_avg[i] < min_25[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_25_a_avg[i] : min_25[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_25[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_25_a_avg[i];

                threshold_30_sk[jj][i] = threshold_30_a_avg[i] > threshold_30[0] ? 65536 : threshold_30_a_avg[i] < threshold_30[1] ? 1 : 0;
                max_30[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_30_a_avg[i] > max_30[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_30_a_avg[i] : max_30[(lwe_sk[i].ConvertToInt()+1) % 65537];
                min_30[(lwe_sk[i].ConvertToInt()+1) % 65537] = threshold_30_a_avg[i] < min_30[(lwe_sk[i].ConvertToInt()+1) % 65537] ? threshold_30_a_avg[i] : min_30[(lwe_sk[i].ConvertToInt()+1) % 65537];
                tot_30[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_30_a_avg[i];
                
                sk_5_cnt[jj] += (threshold_5_sk[jj][i] == lwe_sk[i]);
                sk_10_cnt[jj] += (threshold_10_sk[jj][i] == lwe_sk[i]);
                sk_15_cnt[jj] += (threshold_15_sk[jj][i] == lwe_sk[i]);
                sk_20_cnt[jj] += (threshold_20_sk[jj][i] == lwe_sk[i]);
                sk_25_cnt[jj] += (threshold_25_sk[jj][i] == lwe_sk[i]);
                sk_30_cnt[jj] += (threshold_30_sk[jj][i] == lwe_sk[i]);
                // cout << lwe_sk[i] << " -- " << 
                //     threshold_5_sk[jj][i] << " , " << 
                //     threshold_10_sk[jj][i] << " , " << 
                //     threshold_15_sk[jj][i] << endl;
            }
            max_success[0] = sk_5_cnt[jj] > max_success[0] ? sk_5_cnt[jj] : max_success[0];
            max_success[1] = sk_10_cnt[jj] > max_success[1] ? sk_10_cnt[jj] : max_success[1];
            max_success[2] = sk_15_cnt[jj] > max_success[2] ? sk_15_cnt[jj] : max_success[2];
            max_success[3] = sk_20_cnt[jj] > max_success[3] ? sk_20_cnt[jj] : max_success[3];
            max_success[4] = sk_25_cnt[jj] > max_success[4] ? sk_25_cnt[jj] : max_success[4];
            max_success[5] = sk_30_cnt[jj] > max_success[5] ? sk_30_cnt[jj] : max_success[5];

            avg_success[0] += sk_5_cnt[jj];
            avg_success[1] += sk_10_cnt[jj];
            avg_success[2] += sk_15_cnt[jj];
            avg_success[3] += sk_20_cnt[jj];
            avg_success[4] += sk_25_cnt[jj];
            avg_success[5] += sk_30_cnt[jj];

            min_success[0] = sk_5_cnt[jj] < min_success[0] ? sk_5_cnt[jj] : min_success[0];
            min_success[1] = sk_10_cnt[jj] < min_success[1] ? sk_10_cnt[jj] : min_success[1];
            min_success[2] = sk_15_cnt[jj] < min_success[2] ? sk_15_cnt[jj] : min_success[2];
            min_success[3] = sk_20_cnt[jj] < min_success[3] ? sk_20_cnt[jj] : min_success[3];
            min_success[4] = sk_25_cnt[jj] < min_success[4] ? sk_25_cnt[jj] : min_success[4];
            min_success[5] = sk_30_cnt[jj] < min_success[5] ? sk_30_cnt[jj] : min_success[5];

            // cout << threshold_5_a_avg << endl << threshold_10_a_avg << endl << threshold_15_a_avg << endl;
            cout << max_5 << "            " << max_10 << "            " << max_15 << "            " << max_20 << "            " << max_25 << "            " << max_30 << endl;
            cout << min_5 << "            " << min_10 << "            " << min_15 << "            " << min_20 << "            " << min_25 << "            " << min_30 << endl;
            cout << "Average: \n";
            for (int i = 0; i < 3; i++) {
                cout << tot_5[i] / a_tot[i] << " " << tot_10[i] / a_tot[i] << " " << tot_15[i] / a_tot[i] << " " << tot_20[i] / a_tot[i] << " " << tot_25[i] / a_tot[i] << " " << tot_30[i] / a_tot[i] <<  endl;
            }
            cout << endl;
        }   


        // datafile << "LWE key: \n" << lwe_sk << "\n";
        int same_wrong_5 = 0;
        int same_wrong_10 = 0;
        int same_wrong_15 = 0;
        int same_wrong_20 = 0;
        int same_wrong_25 = 0;
        int same_wrong_30 = 0;
        int all_wrong_5 = 0;
        int all_wrong_10 = 0;
        int all_wrong_15 = 0;
        int all_wrong_20 = 0;
        int all_wrong_25 = 0;
        int all_wrong_30 = 0;

        for (int i = 0; i < n; i++) {
            bool same = true;
            bool all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_5_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_5_sk[j-1][i] == threshold_5_sk[j][i]);
            }
            all_wrong_5 += (int) all_wrong;
            same_wrong_5 += (int) (all_wrong && same);
            
            same = true;
            all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_10_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_10_sk[j-1][i] == threshold_10_sk[j][i]);
            }
            all_wrong_10 += (int) all_wrong;
            same_wrong_10 += (int) (all_wrong && same);

            same = true;
            all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_15_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_15_sk[j-1][i] == threshold_15_sk[j][i]);
            }
            all_wrong_15 += (int) all_wrong;
            same_wrong_15 += (int) (all_wrong && same);

            same = true;
            all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_20_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_20_sk[j-1][i] == threshold_20_sk[j][i]);
            }
            all_wrong_20 += (int) all_wrong;
            same_wrong_20 += (int) (all_wrong && same);

            same = true;
            all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_25_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_25_sk[j-1][i] == threshold_25_sk[j][i]);
            }
            all_wrong_25 += (int) all_wrong;
            same_wrong_25 += (int) (all_wrong && same);

            same = true;
            all_wrong = true; 
            for (int j = 0; j < group_cnt; j++){
                all_wrong = (all_wrong && threshold_30_sk[j][i] != lwe_sk[i]);
                if (j != 0) same = (same && threshold_30_sk[j-1][i] == threshold_30_sk[j][i]);
            }
            all_wrong_30 += (int) all_wrong;
            same_wrong_30 += (int) (all_wrong && same);
        }
        // datafile << "success rate for threshold 2: " << sk_5_cnt << "\n" << "success rate for threshold 3: " << sk_10_cnt << "\n" << "success rate for threshold 5: " << sk_15_cnt << "\n";
        // datafile << "for threshold 2,3,5 count for all wrong: " << all_wrong_5 << ", " << all_wrong_10 << ", " << all_wrong_15 << "\n";
        // datafile << "for threshold 2,3,5 count for all same wrong: " << same_wrong_5 << ", " << same_wrong_10 << ", " << same_wrong_15 << "\n\n";

        max_same_wrong[0] = same_wrong_5 > max_same_wrong[0] ? same_wrong_5 : max_same_wrong[0];
        max_same_wrong[1] = same_wrong_10 > max_same_wrong[1] ? same_wrong_10 : max_same_wrong[1];
        max_same_wrong[2] = same_wrong_15 > max_same_wrong[2] ? same_wrong_15 : max_same_wrong[2];
        max_same_wrong[3] = same_wrong_20 > max_same_wrong[3] ? same_wrong_20 : max_same_wrong[3];
        max_same_wrong[4] = same_wrong_25 > max_same_wrong[4] ? same_wrong_25 : max_same_wrong[4];
        max_same_wrong[5] = same_wrong_30 > max_same_wrong[5] ? same_wrong_30 : max_same_wrong[5];

        avg_same_wrong[0] += same_wrong_5;
        avg_same_wrong[1] += same_wrong_10;
        avg_same_wrong[2] += same_wrong_15;
        avg_same_wrong[3] += same_wrong_20;
        avg_same_wrong[4] += same_wrong_25;
        avg_same_wrong[5] += same_wrong_30;

        max_all_wrong[0] = all_wrong_5 > max_all_wrong[0] ? all_wrong_5 : max_all_wrong[0];
        max_all_wrong[1] = all_wrong_10 > max_all_wrong[1] ? all_wrong_10 : max_all_wrong[1];
        max_all_wrong[2] = all_wrong_15 > max_all_wrong[2] ? all_wrong_15 : max_all_wrong[2];
        max_all_wrong[3] = all_wrong_20 > max_all_wrong[3] ? all_wrong_20 : max_all_wrong[3];
        max_all_wrong[4] = all_wrong_25 > max_all_wrong[4] ? all_wrong_25 : max_all_wrong[4];
        max_all_wrong[5] = all_wrong_30 > max_all_wrong[5] ? all_wrong_30 : max_all_wrong[5];

        avg_all_wrong[0] += all_wrong_5;
        avg_all_wrong[1] += all_wrong_10;
        avg_all_wrong[2] += all_wrong_15;
        avg_all_wrong[3] += all_wrong_20;
        avg_all_wrong[4] += all_wrong_25;
        avg_all_wrong[5] += all_wrong_30;
    }

    for (int i = 0; i < (int) avg_success.size(); i++) {
        avg_success[i] = avg_success[i] / group_cnt/skk_cnt;
        avg_all_wrong[i] = avg_all_wrong[i] / skk_cnt;
        avg_same_wrong[i] = avg_same_wrong[i] / skk_cnt;
    }

    cout << "Avg success: " << avg_success << endl;
    cout << "Max success: " << max_success << endl;
    cout << "Min success: " << min_success << endl;
    cout << "Avg all wrong: " << avg_all_wrong << endl;
    cout << "Max all wrong: " << max_all_wrong << endl;
    cout << "Avg all same wrong: " << avg_same_wrong << endl;
    cout << "Max all same wrong: " << max_same_wrong << endl;


    datafile.close();


}
