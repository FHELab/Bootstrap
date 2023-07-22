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


    ////////////////////////////////////////////// PREPARE (R)LWE PARAMS ///////////////////////////////////////////////
    int ring_dim = poly_modulus_degree_glb;
    int n = 4;
    int p = prime_p;
    int sq_rt = 128;


    EncryptionParameters bfv_params(scheme_type::bgv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 60, 55, 28, 60, 60,
                                                          60, 60, 60, 60, 60,
                                                          50, 60 });

    // auto coeff_modulus = CoeffModulus::Create(ring_dim, { 47, 45, 60, 60, 60,
    //                                                       60, 60, 60, 60, 60,
    //                                                       50, 60 }); // 0 even, 2 odd, for 1,2,3,4
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    bfv_params.set_random_generator(rng);

    SEALContext seal_context(bfv_params, true, sec_level_type::none);
    cout << "primitive root: " << seal_context.first_context_data()->plain_ntt_tables()->get_root() << endl;

    KeyGenerator new_key_keygen(seal_context, n);
    SecretKey new_key = new_key_keygen.secret_key();
    inverse_ntt_negacyclic_harvey(new_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);

    for (int i = 0; i < n; i++) {
        cout <<  new_key.data()[i] << ", ";
    }
    cout << endl << endl;
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        // lwe_sk[i] = (uint64_t) new_key.data()[i] > (uint64_t) p ? p-1 : new_key.data()[i];
        lwe_sk[i] = (uint64_t) new_key.data()[i];
    }
    auto lwe_sk_mod = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        lwe_sk_mod[i] = (uint64_t) new_key.data()[i] > (uint64_t) p ? p-1 : new_key.data()[i];
        // lwe_sk[i] = (uint64_t) new_key.data()[i];
        cout <<  lwe_sk_mod[i] << ", ";
    }
    // cout << endl;
    seal::util::RNSIter new_key_rns(new_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    KSwitchKeys ksk;
    seal::util::ConstPolyIter secret_key_before(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    new_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
    ksk.parms_id() = seal_context.key_parms_id();


    MemoryPoolHandle my_pool = MemoryPoolHandle::New();


    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);
    
    GaloisKeys gal_keys, gal_keys_coeff;
    vector<int> rot_steps = {1};

    keygen.create_galois_keys(rot_steps, gal_keys);
    
    vector<Modulus> coeff_modulus_last = coeff_modulus;
    coeff_modulus_last.erase(coeff_modulus_last.begin() + 3, coeff_modulus_last.end()-1);
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
        i += sq_rt;
    }
    KeyGenerator keygen_last(seal_context_last, sk_last);
    keygen_last.create_galois_keys(rot_steps_coeff, gal_keys_coeff);



    Plaintext pl;
    pl.resize(ring_dim);
    pl.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) ring_dim; i++) {
        pl.data()[i] = i %2 == 0 ? 0 : 4;
    }
    
    Ciphertext coeff;
    encryptor.encrypt(pl, coeff);


    while(seal_context.last_parms_id() != coeff.parms_id()){
        evaluator.mod_switch_to_next_inplace(coeff);
    }


    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    vector<int> msg_r(ring_dim);
    vector<regevCiphertext> lwe_ct_results = extractBGV_noMod(coeff);
    regevDec_BGV_Mod_noMod(msg_r, lwe_ct_results, lwe_sk, lwe_params);
    
    cout << "MSG: \n" ;
    for (int i = 0; i < 4; i++) {
        cout << msg_r[i] << ", ";
    }
    cout << endl;


    cout << "Print extracted LWE: " << endl;
    for (int i = 0; i < 4; i++) {
        cout << "A: ";
        for (int j = 0; j < n; j++) {
            cout << lwe_ct_results[i].a[j].ConvertToInt() << ", ";
        }
        cout << endl;
        cout << "B: " << lwe_ct_results[i].b.ConvertToInt() << endl;
    }

    vector<regevCiphertext> final = extractBGV_test(lwe_ct_results, ring_dim);

    cout << "Print final extracted LWE: " << endl;
    for (int i = 0; i < 4; i++) {
        cout << "A: ";
        for (int j = 0; j < n; j++) {
            cout << final[i].a[j].ConvertToInt() << ", ";
        }
        cout << endl;
        cout << "B: " << final[i].b.ConvertToInt() << endl;
    }


    vector<int> msg_rr(ring_dim);
    regevDec_BGV_Mod3(msg_rr, final, lwe_sk_mod, lwe_params);


    cout << "MSG: \n" ;
    for (int i = 0; i < 4; i++) {
        cout << msg_rr[i] << ", ";
    }
    cout << endl;









    // ///////////// TEST NEW RANGE CHECK /////////////////////
    // Ciphertext output;

    // map<int, bool> modDownIndices_1 = {{4, false}, {12, false}};
    // map<int, bool> modDownIndices_2 = {{4, false}, {16, false}};
    
    // chrono::high_resolution_clock::time_point time_start, time_end;
    // time_start = chrono::high_resolution_clock::now();
    // Bootstrap_FastRangeCheck_Condition(bfv_secret_key, output, c, poly_modulus_degree_glb, relin_keys, seal_context, fastRangeCheckIndices_63_bigPrime,
    //                                    8, 16, modDownIndices_1, modDownIndices_2, 256, 256, modDownIndices_1, modDownIndices_1);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // cout << decryptor.invariant_noise_budget(output) << " bits\n";

    // decryptor.decrypt(output, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "MSG ----------------\n" << msg << endl;






    // ///////////////// TEST BIG RING DIM SLOT TO COEFF ////////////////////


    // vector<Ciphertext> ct_sqrt_list(2*sq_ct);

    // for (int i = 0; i < 9; i++) {
    //     evaluator.mod_switch_to_next_inplace(c);
    // }
    // cout << "... prepared bfv input ciphertext nearly out of noise budget ...\n";
    // cout << decryptor.invariant_noise_budget(c) << endl;
    // Ciphertext bfv_input_copy(c);

    // Evaluator eval_coeff(seal_context_last);
    // eval_coeff.rotate_columns_inplace(c, gal_keys_coeff);
    // for (int i = 0; i < sq_ct; i++) {
    //     eval_coeff.rotate_rows(c, sq_rt * i, gal_keys_coeff, ct_sqrt_list[i]);
    //     eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i]);
    //     eval_coeff.rotate_rows(bfv_input_copy, sq_rt * i, gal_keys_coeff, ct_sqrt_list[i+sq_ct]);
    //     eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    // }

    // vector<Plaintext> U_plain_list(ring_dim);
    // for (int iter = 0; iter < sq_rt; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         vector<uint64_t> U_tmp = readUtemp(j*sq_rt + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }


    // cout << "... prepared rotated bfv input ciphertext ...\n";

    // decryptor.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "Decrypted should be: " << msg << endl;

    // Ciphertext coeff = slotToCoeff_bigRingDim(seal_context, seal_context_last, ct_sqrt_list, U_plain_list, gal_keys_coeff, sq_rt, ring_dim);
    // // Ciphertext coeff = slotToCoeff_WOPrepreocess_bigRingDim(seal_context, seal_context_last, ct_sqrt_list, gal_keys_coeff, sq_rt, ring_dim, p);

    // cout << decryptor.invariant_noise_budget(coeff) << endl;

    // decryptor.decrypt(coeff, pl);
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << pl.data()[i] << " ";
    // }
    // cout << endl;









    // Evaluator eval_gal(seal_context_last);

    // Ciphertext c_mod(c);
    // evaluator.mod_switch_to_next_inplace(c_mod);

    // cout << "before mod...\n";
    // while (seal_context.last_parms_id() != c.parms_id()) {
    //     // evaluator.mod_switch_to_next_inplace(c_mod);
    //     evaluator.mod_switch_to_next_inplace(c);
    // }

    // int sq_sk = sqrt(n);
    // Ciphertext c_column;
    // vector<Ciphertext> lwe_sk_sqrt_list(sq_sk);
    // cout << "1\n";
    // eval_gal.rotate_columns(c, gal_keys_coeff, c_column);
    // for (int i = 0; i < sq_sk; i++) {
    //     cout << i << endl;
    //     eval_gal.rotate_rows(c, sq_sk * i, gal_keys_coeff, lwe_sk_sqrt_list[i]);
    // }


    // Plaintext pl_1;

    // Ciphertext c2;
    // chrono::high_resolution_clock::time_point time_start, time_end;
    // time_start = chrono::high_resolution_clock::now();
    // for (int i = 0; i < ring_dim; i++) {
    //     evaluator.multiply_plain(c, pl, c2);
    // }
    // time_end = chrono::high_resolution_clock::now();
    // cout << "Pl Multi: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // pl_1.resize(ring_dim);
    // pl_1.parms_id() = parms_id_zero;
    // cout << "?\n";
    // pl_1.data()[0] = 2;
    // for (int i = 1; i < ring_dim; i++) {
    //     pl_1.data()[i] = 0;
    // }

    // cout << "Before multi..\n";
    // evaluator.multiply_plain_inplace(c, pl_1);

    // decryptor.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "Decrypted : " << msg << endl;


    // Ciphertext c1;
    // encryptor.encrypt(pl_1, c1);
    // time_start = chrono::high_resolution_clock::now();
    // for (int i = 0; i< 10; i++) {
    //     evaluator.multiply_inplace(c, c1);
    // }
    // time_end = chrono::high_resolution_clock::now();
    // cout << "Ct Multi: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;




    // TEST KEY SWITCH.....


    // auto coeff_modulus_switch = CoeffModulus::Create(ring_dim, { 28, 60 });
    // EncryptionParameters parms_switch = bfv_params;
    // parms_switch.set_coeff_modulus(coeff_modulus_switch);
    // SEALContext seal_context_switch = SEALContext(parms_switch, true, sec_level_type::none);


    // KeyGenerator new_key_keygen(seal_context_switch, n);
    // SecretKey new_key = new_key_keygen.secret_key();
    // KSwitchKeys ksk;
    // seal::util::ConstPolyIter secret_key_before(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    // new_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
    // ksk.parms_id() = seal_context.key_parms_id();

    // Evaluator eval_switch(seal_context_switch);







    // decryptor.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "Original Decrypt: " << msg << endl;



    // cout << "New SK: " << endl;
    // for (int i = 0; i < 10; i++) {
    //     // cout << new_key.data()[i] << " --> ";
    //     // new_key.data()[i] = (uint64_t) ((float(new_key.data()[i])) * float(268369921) / float(1152921504581419009));
    //     cout << new_key.data()[i] << "  ";
    // }
    // cout << endl;

    // while(seal_context.last_parms_id() != c.parms_id()){
    //     evaluator.mod_switch_to_next_inplace(c);
    // }


    // // 36028797017456641
    // // 1152921504606584833 / 32990138759887301

    // modDownToPrime(c, ring_dim, 1152921504581419009, 268369921);
    // // for (int i = 0; i < ring_dim; i++) {
    // //     cout << c.data(0)[i] << "  ";
    // // }
    // // cout << endl;

    // Ciphertext copy_coeff = c;
    // auto ct_in_iter = util::iter(copy_coeff);
    // ct_in_iter += c.size() - 1;
    // seal::util::set_zero_poly(ring_dim, 1, c.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    // c.parms_id_ = seal_context_switch.last_parms_id();

    // eval_switch.switch_key_inplace(c, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);


    // for (int i = 0; i < ring_dim; i++) {
    //     cout << c.data(0)[i] << "  ";
    // }
    // cout << endl;
    
    // Decryptor decryptor_new(seal_context_switch, new_key);


    // decryptor_new.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "New Decrypt: " << msg << endl;




    // TEST MODDOWN FUNC.....



    // vector<Ciphertext> output(1023);
    // // map<int, bool> modDownIndices = {{4, false}, {16, false}, {64, false}, {256, false}};

    // map<int, bool> modDownIndices = {{2, false}, {8, false}, {32, false}, {64, false}, {128, false}, {512, false}};
    // // chrono::high_resolution_clock::time_point time_start, time_end;
    // // time_start = chrono::high_resolution_clock::now();
    // calUptoDegreeK_bigPrime(output, c, 1023, relin_keys, seal_context, modDownIndices);


    // cout << "Decryted check: \n";
    // vector<uint64_t> v(ring_dim);
    // for (int i = 0; i < output.size(); i++) {
    //     decryptor.decrypt(output[i], pl);
    //     batch_encoder.decode(pl, v);
    //     cout << v[1] << ",";
    // }
    // cout << endl;





    // TEST SLOTTOCOEFF STUFF.....

    // Plaintext pl;
    // Ciphertext c;
    // batch_encoder.encode(msg, pl);

    // for (int i = 0; i < ring_dim; i++) {
    //     cout << pl[i] << " ";
    // }
    // cout << endl;
    // encryptor.encrypt(pl, c);


    // Ciphertext c_copy(c);

    // int sq_ct = sqrt(ring_dim/2);
    // vector<Ciphertext> ct_sqrt_list(2*sq_ct);

    // evaluator.rotate_columns_inplace(c_copy, gal_keys);
    // for (int i = 0; i < sq_ct; i++) {
    //     evaluator.rotate_rows(c, sq_ct * i, gal_keys, ct_sqrt_list[i]);
    //     evaluator.transform_to_ntt_inplace(ct_sqrt_list[i]);
    //     evaluator.rotate_rows(c_copy, sq_ct * i, gal_keys, ct_sqrt_list[i+sq_ct]);
    //     evaluator.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    // }

    // evaluator.rotate_rows_inplace(c, 1, gal_keys);

    // decryptor.decrypt(c, pl);
    // batch_encoder.decode(pl, msg);
    // cout << "Decode: " << msg << endl;


    // for (int i = 0; i < ct_sqrt_list.size(); i++) {
    //     Plaintext pp;
    //     vector<uint64_t> v(ring_dim);

    //     evaluator.transform_from_ntt_inplace(ct_sqrt_list[i]);

    //     decryptor.decrypt(ct_sqrt_list[i], pp);
    //     batch_encoder.decode(pp, v);
    //     cout << v << endl;

    // }

    // vector<Plaintext> U_plain_list(ring_dim);
    // vector<uint64_t> U_tmp;
    // for (int iter = 0; iter < sq_ct; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         U_tmp = readUtemp(j*sq_ct + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }


    // Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, ct_sqrt_list, gal_keys, ring_dim);

    // // evaluator.rotate_columns_inplace(c, gal_keys);

    // decryptor.decrypt(coeff, pl);
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << pl[i] << " ";
    // }
    // cout << endl;

}
