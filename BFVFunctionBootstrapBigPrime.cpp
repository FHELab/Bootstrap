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
    int n = 1024;
    BootstrapParam bootstrap_param = BootstrapParam(prime_p, 192, 4096, 256*3, 1024);
    int p = bootstrap_param.ciphertextSpacePrime;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);

    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 60, 60, 28, 60, 60, 60, 60,
                                                          60, 60, 60, 60, 60, 60, 60,
                                                          50, 60 });
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

    KeyGenerator lwe_key_keygen(seal_context, n);
    SecretKey lwe_key = lwe_key_keygen.secret_key();

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    // generate a key switching key based on key_before and secret_key
    KSwitchKeys ksk_to_lwe;
    seal::util::ConstPolyIter secret_key_bfv(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    lwe_key_keygen.generate_kswitch_keys(secret_key_bfv, 1, static_cast<KSwitchKeys &>(ksk_to_lwe), false); // used to switch from secret_key_bfv to secret_key_lwe
    ksk_to_lwe.parms_id() = seal_context.key_parms_id();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    GaloisKeys gal_keys, gal_keys_coeff;
    vector<int> rot_steps = {1};
    for (int i = 0; i < n;) {
        rot_steps.push_back(i);
        i += sqrt(n);
    }
    keygen.create_galois_keys(rot_steps, gal_keys);
    
    vector<Modulus> coeff_modulus_last = coeff_modulus;
    coeff_modulus_last.erase(coeff_modulus_last.begin() + 1, coeff_modulus_last.end()-2);
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

    ////////////////////////////////////////////// PREPARE BFV CIPHERTEXT //////////////////////////////////////////////

    Ciphertext bfv_input;
    vector<uint64_t> input_v(poly_modulus_degree_glb);
    for (int i = 0; i < (int) poly_modulus_degree_glb; i++) {
        input_v[i] = i % 2 == 0 ? 10 : 10001;
    }
    Plaintext pl;
    batch_encoder.encode(input_v, pl);
    encryptor.encrypt(pl, bfv_input);

    cout << "... prepared bfv input ciphertext ...\n";

    ////////////////////////////////////////////// ENCRYPT SK UNDER BFV ////////////////////////////////////////////////

    // one switching key for one lwe_sk
    // Ciphertext lwe_sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);
    inverse_ntt_negacyclic_harvey(lwe_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);
    auto lwe_params = regevParam(n, p, 1.3, ring_dim); 
    auto lwe_sk = regevGenerateSecretKey(lwe_params);
    for (int i = 0; i < n; i++) {
        lwe_sk[i] = (uint64_t) lwe_key.data()[i] > (uint64_t) p ? p-1 : lwe_key.data()[i];
    }

    seal::util::RNSIter new_key_rns(lwe_key.data().data(), ring_dim);
    ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());

    Ciphertext sk_encrypted = encryptLWEskUnderBFV(seal_context, ring_dim, bfv_public_key, bfv_secret_key, lwe_sk, lwe_params);

    cout << "... prepared lwe sk under bfv for extraction ...\n";

    /////////////////////////////////////////////////// BOOTSTRAP //////////////////////////////////////////////////////
    Evaluator evaluator(seal_context);
    vector<uint64_t> q_shift_constant(ring_dim, 0);

    int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(sk_encrypted, sq_sk * i, gal_keys, sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(sk_sqrt_list[i]);
    }

    for (int i = 0; i < 13; i++) {
        evaluator.mod_switch_to_next_inplace(bfv_input);
    }
    cout << "... prepared bfv input ciphertext nearly out of noise budget ...\n";
    Ciphertext bfv_input_copy(bfv_input);

    Evaluator eval_coeff(seal_context_last);
    eval_coeff.rotate_columns_inplace(bfv_input, gal_keys_coeff);
    for (int i = 0; i < sq_ct; i++) {
        eval_coeff.rotate_rows(bfv_input, sq_ct * i, gal_keys_coeff, ct_sqrt_list[i]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i]);
        eval_coeff.rotate_rows(bfv_input_copy, sq_ct * i, gal_keys_coeff, ct_sqrt_list[i+sq_ct]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    }

    cout << "... prepared rotated bfv input ciphertext ...\n";

    // vector<Plaintext> U_plain_list(ring_dim);
    // for (int iter = 0; iter < sq_ct; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         vector<uint64_t> U_tmp = readUtemp(j*sq_ct + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }



    chrono::high_resolution_clock::time_point time_start, time_end;

    time_start = chrono::high_resolution_clock::now();


    // Ciphertext coeff = slotToCoeff(seal_context, seal_context_last, ct_sqrt_list, U_plain_list, gal_keys_coeff, ring_dim);
    Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, seal_context_last, ct_sqrt_list, gal_keys_coeff, ring_dim, p);

    // time_end = chrono::high_resolution_clock::now();

    cout << "slot noise: " << decryptor.invariant_noise_budget(coeff) << endl;






    while(seal_context.last_parms_id() != coeff.parms_id()){
        evaluator.mod_switch_to_next_inplace(coeff);
    }

    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk_to_lwe), 0, my_pool);





    vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff);

    // vector<int> msg(ring_dim);
    // regevDec_Value(msg, lwe_ct_results, lwe_sk, lwe_params, bootstrap_param.errorRange);






    Ciphertext eval_result = evaluateExtractedBFVCiphertext(seal_context, lwe_ct_results, sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, false);

    decryptor.decrypt(eval_result, pl);
    batch_encoder.decode(pl, input_v);
    cout << "Result after eval with lwe key: ---------------------\n" << input_v << endl;


    Ciphertext range_check_res;
    map<int, bool> modDownIndices = {{4, false}, {16, false}};
    Bootstrap_FastRangeCheck_Condition(bfv_secret_key, range_check_res, eval_result, ring_dim, relin_keys, seal_context, fastRangeCheckIndices_63_bigPrime,
                                    bootstrap_param.firstLevelDegree, bootstrap_param.secondLevelDegree, modDownIndices, modDownIndices);

    cout << "final: " << decryptor.invariant_noise_budget(range_check_res) << endl;

    time_end = chrono::high_resolution_clock::now();



    decryptor.decrypt(range_check_res, pl);
    batch_encoder.decode(pl, input_v);
    cout << "Result !!!!! ---------------------\n" << input_v << endl;

    cout << "TOTAL TIME: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;


}