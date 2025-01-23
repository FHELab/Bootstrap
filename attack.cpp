#include "regevEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace seal;
using namespace std;

// vector<double> gen_hint(Ciphertext ct, int index, int n = 1024) {
//     vector<double> hint(n, 0);

//     hint[0] = ct.data(0)[index]; // b

//     for (int i = 0; i < n; i++) {

//     }
// }


// used to bootstrap for BFV ciphertexts, encrypting [0, t-1, r], where t is the prime, and r is the interval (≥ error bound = 128 in our case)
int main() {

    ////////////////////////////////////////////// PREPARE (R)LWE PARAMS ///////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 1024;
    BootstrapParam bootstrap_param = BootstrapParam(65537, 128, 512, 256, 256);
    int p = bootstrap_param.ciphertextSpacePrime;
    // int interval = 128;
    vector<uint64_t> rangeCheckIndices = rangeCheckIndices_bfv;
    // int scalar = bootstrap_param.errorRange/interval;

    // vector<float> threshold_5 = {0.00582, -0.00822};
    vector<float> threshold_5 = {0.007, -0.007};
    
    int group_cnt = 1;
    int threshold_value_input = 58; // threshold: 5

    vector<float> threshold_a(n, 0);
    int a_cnt = 0;

    diff_global.resize(ring_dim);

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);

    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 60, 55, 60, 60,
                                                        60, 60, 60, 60,
                                                        60, 50, 60 });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    prng_seed_type seed_lwe;
    for (auto &i : seed_lwe) {
        i = random_uint64();
    }
    auto rng_lwe = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed_lwe));
    bfv_params.set_random_generator(rng_lwe);

    SEALContext seal_contex_lwe(bfv_params, true, sec_level_type::none);

    KeyGenerator lwe_key_keygen(seal_contex_lwe, n);
    SecretKey lwe_key = lwe_key_keygen.secret_key();


    inverse_ntt_negacyclic_harvey(lwe_key.data().data(), seal_contex_lwe.key_context_data()->small_ntt_tables()[0]);
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = (uint64_t) lwe_key.data()[i] > (uint64_t) p ? p-1 : lwe_key.data()[i];
        if ((uint64_t) lwe_key.data()[i] > (uint64_t) p) big_prime_global = (uint64_t) lwe_key.data()[i]+1; 
        // cout << lwe_key.data()[i] << " ";
    }
    // cout << "------------------------------\n" << lwe_sk << "\n------------------------------\n";

    seal::util::RNSIter new_key_rns(lwe_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_contex_lwe.key_context_data()->small_ntt_tables());

    for (int gg = 0; gg < group_cnt; gg++) {

        prng_seed_type seed;
        for (auto &i : seed) {
            i = random_uint64();
        }
        auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
        bfv_params.set_random_generator(rng);

        SEALContext seal_context(bfv_params, true, sec_level_type::none);

        KeyGenerator keygen(seal_context);
        SecretKey bfv_secret_key = keygen.secret_key();

        // generate a key switching key based on key_before and secret_key
        KSwitchKeys ksk_to_lwe, ksk_to_bfv;
        seal::util::ConstPolyIter secret_key_bfv(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());
        seal::util::ConstPolyIter secret_key_lwe(lwe_key.data().data(), ring_dim, coeff_modulus.size());

        lwe_key_keygen.generate_kswitch_keys(secret_key_bfv, 1, static_cast<KSwitchKeys &>(ksk_to_lwe), false); // used to switch from secret_key_bfv to secret_key_lwe
        keygen.generate_kswitch_keys(secret_key_lwe, 1, static_cast<KSwitchKeys &>(ksk_to_bfv), false); // used to switch from secret_key_lwe to secret_key_bfv
        ksk_to_lwe.parms_id() = seal_context.key_parms_id();
        ksk_to_bfv.parms_id() = seal_context.key_parms_id();

        PublicKey bfv_public_key;
        keygen.create_public_key(bfv_public_key);

        RelinKeys relin_keys;
        keygen.create_relin_keys(relin_keys);

        Encryptor encryptor(seal_context, bfv_public_key);
        BatchEncoder batch_encoder(seal_context);
        Decryptor decryptor(seal_context, bfv_secret_key);
        Decryptor decryptor_lwe(seal_context, lwe_key);
        GaloisKeys gal_keys, gal_keys_coeff;
        vector<int> rot_steps = {1};
        for (int i = 0; i < n;) {
            rot_steps.push_back(i);
            i += sqrt(n);
        }
        keygen.create_galois_keys(rot_steps, gal_keys);
        
        vector<Modulus> coeff_modulus_last = coeff_modulus;
        coeff_modulus_last.erase(coeff_modulus_last.begin() + 1, coeff_modulus_last.end()-1);
        EncryptionParameters parms_last = bfv_params;
        parms_last.set_coeff_modulus(coeff_modulus_last);
        SEALContext seal_context_last = SEALContext(parms_last, true, sec_level_type::none);

        SecretKey sk_last;
        sk_last.data().resize(coeff_modulus_last.size() * ring_dim);
        sk_last.parms_id() = seal_context_last.key_parms_id();
        util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_last.size() - 1, sk_last.data().data());
        util::set_poly(
            bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
            sk_last.data().data() + ring_dim * (coeff_modulus_last.size() - 1));

        vector<int> rot_steps_coeff = {1};
        for (int i = 0; i < ring_dim/2;) {
            if (find(rot_steps_coeff.begin(), rot_steps_coeff.end(), i) == rot_steps_coeff.end()) {
                rot_steps_coeff.push_back(i);
            }
            i += sqrt(ring_dim/2);
        }
        KeyGenerator keygen_last(seal_context_last, sk_last);
        keygen_last.create_galois_keys(rot_steps_coeff, gal_keys_coeff);

        MemoryPoolHandle my_pool = MemoryPoolHandle::New();

        // cout << "After param gen.\n";

        ////////////////////////////////////////////// PREPARE BFV CIPHERTEXT //////////////////////////////////////////////

        // Ciphertext bfv_input;
        vector<uint64_t> input_v(poly_modulus_degree_glb);
        // for (int i = 0; i < (int) poly_modulus_degree_glb; i++) {
        //     // input_v[i] = (i % 512) * interval;
        //     input_v[i] = 54;
        // }
        // Plaintext pl;
        // batch_encoder.encode(input_v, pl);
        // encryptor.encrypt(pl, bfv_input);

        Ciphertext bfv_input;
        Plaintext pl;
        pl.resize(poly_modulus_degree_glb);
        pl.parms_id() = parms_id_zero;
        // vector<uint64_t> input_v(poly_modulus_degree_glb);
        for (int i = 0; i < (int) poly_modulus_degree_glb; i++) {
            pl.data()[i] = threshold_value_input;
        }
        // batch_xencoder.encode(input_v, pl);
        encryptor.encrypt(pl, bfv_input);

        ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

        // one switching key for one lwe_sk
        // Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);

        Ciphertext sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);


        /////////////////////////////////////////////////// BOOTSTRAP //////////////////////////////////////////////////////
        Evaluator evaluator(seal_context);
        vector<uint64_t> q_shift_constant(ring_dim, 0);

        int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
        vector<Ciphertext> sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

        for (int i = 0; i < sq_sk; i++) {
            evaluator.rotate_rows(sk_encrypted, sq_sk * i, gal_keys, sk_sqrt_list[i]);
            evaluator.transform_to_ntt_inplace(sk_sqrt_list[i]);
        }

        for (int i = 0; i < 9; i++) {
            evaluator.mod_switch_to_next_inplace(bfv_input);
        }


        while(seal_context.last_parms_id() != bfv_input.parms_id()){
            evaluator.mod_switch_to_next_inplace(bfv_input);
        }

        // decryptor.decrypt(coeff, pl);
        // for (int i = 0; i < (int) poly_modulus_degree_glb; i++) {
        //     cout << pl.data()[i] << " ";
        // }
        // cout << "hahaha\n" << endl;


        Ciphertext copy_coeff = bfv_input;
        auto ct_in_iter = util::iter(copy_coeff);
        ct_in_iter += bfv_input.size() - 1;
        seal::util::set_zero_poly(ring_dim, 1, bfv_input.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level
        evaluator.switch_key_inplace(bfv_input, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk_to_lwe), 0, my_pool);


        // decryptor_lwe.decrypt(coeff, pl);
        // for (int i = 0; i < (int) ring_dim; i++) {
        // cout << pl.data()[i] << " ";
        // }
        // cout << endl << endl;
        vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(bfv_input, ring_dim, n);


        vector<int> ttmsg;
        regevDec(ttmsg, lwe_ct_results, lwe_sk, lwe_params);

        // cout << "heihei\n" << ttmsg << endl;


        Ciphertext eval_result = evaluateExtractedBFVCiphertext(seal_context, lwe_ct_results, sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, false);

        // decryptor.decrypt(eval_result, pl);
        // batch_encoder.decode(pl, input_v);
        // cout << "Result after eval with lwe key: ---------------------\n" << input_v << endl;


        Ciphertext range_check_res;
        /* for bootstrap function evaluation, use rangeCheckIndices_bfv for identity mapping, and rangeCheckIndices_bfv_square for square mapping */
        Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, eval_result, rangeCheckIndices, p, ring_dim,
                                                relin_keys, seal_context, bfv_secret_key, 0, false, false,
                                                bootstrap_param.firstLevelDegree, bootstrap_param.secondLevelDegree);


        decryptor.decrypt(range_check_res, pl);
        batch_encoder.decode(pl, input_v);
        // cout << "Result !!!!! ---------------------\n" << input_v << endl;
        

        for (int mi = 0; mi < (int) input_v.size(); mi++) {
            if (input_v[mi] == 1) {
                a_cnt+=1;
                for (int ai = 0; ai < n; ai++) {
                    threshold_a[ai] += diff_global[mi][ai];
                }
            }
        }
    }

    cout << "Collect " << a_cnt << " sampled.\n";

    vector<int> threshold_sk(n, 0);
    int sk_cnt = 0;
    vector<float> tot(3, 0);
    vector<int> a_tot(3, 0);

    for (int i = 0; i < n; i++) {
        threshold_a[i] = threshold_a[i] / a_cnt;
    
        threshold_sk[i] = threshold_a[i] > threshold_5[0] ? 65536 : threshold_a[i] < threshold_5[1] ? 1 : 0;
        
        sk_cnt += (threshold_sk[i] == lwe_sk[i]);
        tot[(lwe_sk[i].ConvertToInt()+1) % 65537] += threshold_a[i];
        a_tot[(lwe_sk[i].ConvertToInt()+1) % 65537] += 1;
        // cout << lwe_sk[i] << " -- " << 
        //     threshold_5_sk[jj][i] << " , " << 
        //     threshold_10_sk[jj][i] << " , " << 
        //     threshold_15_sk[jj][i] << endl;
    }

    // cout << "Threshold values: " << threshold_a << endl;

    // cout << lwe_sk << endl << threshold_sk << endl;
    cout << "success rate: " << sk_cnt << endl;

    for (int i = 0; i < 3; i++) {
        cout << tot[i] / a_tot[i] << endl;
    }



}
