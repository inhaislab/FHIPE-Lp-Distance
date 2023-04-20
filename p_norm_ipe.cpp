//
// Created by DongHyeon Yu on 22. 11. 28..
//

#include "p_norm_ipe.h"

#include <cmath>

#include "ffhipe.h"
#include "ternary_encoder.hpp"
#include "dummy_encoder.hpp"
#include "msk.h"
#include "util.h"
#include "fhipe_manager.h"
#include "p_norm_encoder.hpp"


/**
 * @param n the dimension of vector X and Y
 * @param p degree. exactly same with "p" when we say "p-power mean error" or "Lp norm".
 */
ffhipe::Msk setup(ULONG n, ULONG p)
{
    ULONG encoded_dimension = (p - 1) * n + 2;
    return ffhipe::FHIPEManager::setup(encoded_dimension, ffhipe::CPubParam::getInstance());
}


/**
 * @param ct an encoded normal state vector with encrypted.
 * @param sk an encoded current state vector with encrypted.
 * @param m max range of each element of vectors.
 * @return Lp norm between normal state vector and current state vector.
 */
MSG decrypt(ffhipe::Ct ct, ffhipe::Sk sk, ULONG m, ULONG p)
{
    ffhipe::CPubParam* param = ffhipe::CPubParam::getInstance();

    ULONG mp = 1;
    for(size_t i = 0; i < p; i++) mp *= m;
    ULONG n = (sk->dim()-2)/(p-1);
    MSG p_powered_inner_prod_value = param->decrypt(sk, ct, mp * n);

    return p_powered_inner_prod_value;
}

void p_norm_ipe::decrypt_value_test(ULONG n, ULONG m, ULONG p) {

    auto encoder = new ffhipe::PNormEncoder<long>(p);

    ffhipe::Msk msk = setup(n, p);

    srand(time(nullptr));

    std::vector<long> _x(n);
    std::vector<long> _y(n);

    MSG solution;
    LONG temp_x, temp_y, temp_diff;

    solution = 0;
    for (ULONG i = 0; i < n; i++) {

        temp_x = rand() % m;
        temp_y = rand() % m;
        _x[i] = temp_x;
        _y[i] = temp_y;

        temp_diff = temp_x - temp_y;
        solution += std::pow(temp_diff, p);
    }

    std::cout << "_x : ";
    for (auto i = 0; i < n; i++) {
        std::cout << _x[i] << ", ";
    }
    std::cout << std::endl;

    std::cout << "_y : ";
    for (auto i = 0; i < n; i++) {
        std::cout << _y[i] << ", ";
    }
    std::cout << std::endl;

    VECTOR<MSG> x, y;

    auto sk = msk->encrypt(encoder->encodeX(_x, x));
    auto ct = msk->keygen(encoder->encodeY(_y, y));

    auto k1 = sk->getM1();
    auto k2 = sk->getM2();
    auto c1 = ct->getM1();
    auto c2 = ct->getM2();

    // It would be good to normalize K1, K2, C1 and C2 when required to transfer data
    k1.normalize();
    for (size_t i = 0; i < sk->dim(); i++) {
        k2[i].normalize();
    }
    c1.normalize();
    for (size_t i = 0; i < ct->dim(); i++) {
        c2[i].normalize();
    }

    MSG z = decrypt(sk, ct, m, p);

    std::cout << "z(" << z << ") is " << solution << "? " << ((z == solution) ? "Yes" : "No") << std::endl;
    std::cout << std::endl;

    auto g1 = msk->getG1();
    auto g2 = msk->getG2();

    std::cout << "g1    : " << g1.x << " " << g1.y << std::endl;
    std::cout << "g2    : " << g2.x << " " << g2.y << std::endl;
    std::cout << std::endl;

    std::cout << "k1    : " << k1.x << " " << k1.y << std::endl;
    for (size_t i = 0; i< sk->dim(); i++) {
        std::cout << "k2(" << i << ") : "<< k2[i].x << " " << k2[i].y << std::endl;
    }
    std::cout << std::endl;

    std::cout << "c1    : " << c1.x << " " << c1.y << std::endl;
    for (size_t i = 0; i< sk->dim(); i++) {
        std::cout << "c2(" << i << ") : "<< c2[i].x << " " << c2[i].y << std::endl;
    }
    std::cout << std::endl;

    delete sk;
    delete ct;

    delete encoder;
    delete msk;
}

#include <fstream>
#include <sstream>

void get_msk_from_file(unsigned long n, unsigned long m, unsigned long &p, ffhipe::Msk &msk) {
    std::stringstream filename;
    filename << "N=" << n << ", M=" << m << ", P=" << p << ".txt";
    std::ifstream ifs(filename.str(), std::ios::in);
    msk = ffhipe::CMsk::load(ifs, p);
}

void p_norm_ipe::encrypt(ULONG n, ULONG m, ULONG p) {

    auto encoder = new ffhipe::PNormEncoder<long>(p);

    srand(time(nullptr));

    std::vector<long> _x(n);
    for (ULONG i = 0; i < n; i++)
        _x[i] = rand() % m;

    VECTOR<MSG> x;

    ffhipe::Msk msk;
    get_msk_from_file(n, m, p, msk);
    delete msk->encrypt(encoder->encodeX(_x, x));
    delete msk;
    delete encoder;
}


void p_norm_ipe::keygen(ULONG n, ULONG m, ULONG p) {

    auto encoder = new ffhipe::PNormEncoder<long>(p);

    srand(time(nullptr));

    std::vector<long> _y(n);
    for (ULONG i = 0; i < n; i++)
        _y[i] = rand() % m;

    VECTOR<MSG> y;

    ffhipe::Msk msk;
    get_msk_from_file(n, m, p, msk);
    delete msk->keygen(encoder->encodeY(_y, y));
    delete msk;
    delete encoder;
}

enum CheckPoint {
    IPE_Setup,
    IPE_ENCODE_X,
    IPE_ENCODE_Y,
    IPE_Keygen,
    IPE_Encrypt,
    IPE_Decrypt,
    NUM_OF_CHECK_POINT
};

void test_p_norm_ipe(ULONG n, ULONG m, ULONG p, ffhipe::tclock::TInfo ti[]) {

    auto * encoder = new ffhipe::PNormEncoder<long>(p);

    ffhipe::tclock::start(ti[IPE_Setup]);
    ffhipe::Msk msk = setup(n, p);
    ffhipe::tclock::end(ti[IPE_Setup]);

    srand(time(nullptr));

    std::vector<long> _x(n);
    std::vector<long> _y(n);

    MSG solution;
    LONG temp_x, temp_y, temp_diff;

    solution = 0;
    for (ULONG i = 0; i < n; i++) {

        temp_x = rand() % m;
        temp_y = rand() % m;
        _x[i] = temp_x;
        _y[i] = temp_y;

        temp_diff = temp_x - temp_y;
        solution += std::pow(temp_diff, p);
    }

    VECTOR<MSG> x, y;

    ffhipe::tclock::start(ti[IPE_ENCODE_X]);
    VECTOR<MSG> encoded_X = encoder->encodeX(_x, x);
    ffhipe::tclock::end(ti[IPE_ENCODE_X]);

    ffhipe::tclock::start(ti[IPE_ENCODE_Y]);
    VECTOR<MSG> encoded_Y = encoder->encodeY(_y, y);
    ffhipe::tclock::end(ti[IPE_ENCODE_Y]);

    ffhipe::tclock::start(ti[IPE_Encrypt]);
    auto sk = msk->encrypt(encoded_X);
    ffhipe::tclock::end(ti[IPE_Encrypt]);

    ffhipe::tclock::start(ti[IPE_Keygen]);
    auto ct = msk->keygen(encoded_Y);
    ffhipe::tclock::end(ti[IPE_Keygen]);

    ffhipe::tclock::start(ti[IPE_Decrypt]);
    MSG z = decrypt(sk, ct, m, p);
    ffhipe::tclock::end(ti[IPE_Decrypt]);

    delete sk;
    delete ct;

    delete encoder;
    delete msk;
}

void p_norm_ipe::decrypt_performance_test(const std::initializer_list<ULONG>& n_arr, ULONG m,const std::initializer_list<ULONG>& p_arr, ULONG times) {

    std::cout << "=== P-norm IPE Performance Measurement ===" << std::endl;
    ffhipe::tclock::TInfo ti[NUM_OF_CHECK_POINT]{};

    std::cout << "+ Details for " << times << "times" << std::endl;
    for (ULONG n : n_arr) {
        for (ULONG p : p_arr)
        {
            std::cout << std::endl << "N=" << n << ", M=" << m  << ", P=" << p << std::endl;
            for (int i = 0; i < times; i++) {
                ffhipe::tclock::init(ti, NUM_OF_CHECK_POINT);
                test_p_norm_ipe(n, m, p, ti);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_Setup]);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_ENCODE_X]);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_ENCODE_Y]);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_Encrypt]);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_Keygen]);
                ffhipe::tclock::outAccTime(std::cout, ti[IPE_Decrypt], true);
            }
        }
    }
}

void get_latent_vector_from_file(int n, std::vector<long>& _dst)
{
    std::stringstream filename;
    filename << "data/latent_vectors/normal/quatized_latent_vectors_" << n*4 << "_" << n*2 << "_" << n ;

    static std::iostream::pos_type _latent_vector_file_offset = 0;
    std::ifstream ofs(filename.str(), std::ios::in);
    if(!ofs.is_open()) {
        std::cerr << "[p_norm_ipe::get_latent_vector_from_file] file isn't opened: " << filename.str() << "\n";
        exit(1);
    }
    ofs.seekg(_latent_vector_file_offset, std::ios::beg);

    //load a latent vector from file.
    std::string line;
    std::getline(ofs, line);
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ' ')) {
        if(!cell.empty())
            _dst.push_back(std::stol(cell));
    }

    if (_dst.size() != n)
    {
        std::cerr << "[p_norm_ipe::get_latent_vector_from_file] no matching vector size: " << _dst.size() << ".\n";
        exit(1);
    }

    _latent_vector_file_offset = ofs.tellg();
    ofs.close();
}

void test_encrypt_keygen_only(ULONG n, ULONG m, ULONG p, ffhipe::tclock::TInfo ti[])
{
    auto * encoder = new ffhipe::PNormEncoder<long>(p);

    std::vector<long> _x;
    get_latent_vector_from_file(n, _x);
    std::vector<long> _y = _x;

    //load a msk from file.
    ffhipe::Msk msk;
    get_msk_from_file(n, m, p, msk);

    VECTOR<MSG> x, y;

    ffhipe::tclock::start(ti[IPE_Encrypt]);
    auto sk = msk->encrypt(encoder->encodeX(_x, x));
    ffhipe::tclock::end(ti[IPE_Encrypt]);

    ffhipe::tclock::start(ti[IPE_Keygen]);
    auto ct = msk->keygen(encoder->encodeY(_y, y));
    ffhipe::tclock::end(ti[IPE_Keygen]);

    delete sk;
    delete ct;
    delete msk;
}

void p_norm_ipe::encrypt_keygen_performance_test(const std::initializer_list<ULONG>& n_arr, ULONG m, const std::initializer_list<ULONG>& p_arr, ULONG times)
{
    std::cout << "=== P-norm IPE encrypt&keygen Performance Measurement ===" << std::endl;
    std::cout << "+ Details for " << times << "times" << std::endl;

    ffhipe::tclock::TInfo ti[NUM_OF_CHECK_POINT]{};

    for (ULONG n : n_arr) {
        for (ULONG p : p_arr)
        {
            ffhipe::tclock::init(ti, NUM_OF_CHECK_POINT);
            std::cout << std::endl << "N=" << n << ", M=" << m  << ", P=" << p << std::endl;

            for (int i = 0; i < times; i++) {

                test_encrypt_keygen_only(n, m, p, ti);
            }

            ffhipe::tclock::outAverageTime(std::cout, ti[IPE_Encrypt]);
            ffhipe::tclock::outAverageTime(std::cout, ti[IPE_Keygen], true);
        }
    }
}



