// @author: Miroslav Puskaric https://github.com/mpuskaric
//

#include "ciphertext-ser.h"
#include "cryptocontext-ser.h"
#include "key/key-ser.h"
#include "scheme/ckksrns/ckksrns-ser.h"
#include "openfhe.h"
#include <random>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono> 
#include <omp.h> 

using namespace lbcrypto;
using namespace std::chrono;

using CC = CryptoContext<DCRTPoly>;         // crypto context
using CT = Ciphertext<DCRTPoly>;            // ciphertext
using PT = Plaintext;                       // plaintext
using Key_Pair = KeyPair<DCRTPoly>;         // secret/public key par.
using Eval_Key = EvalKey<DCRTPoly>;         // evaluation key (reencryption key)
using Private_Key = PrivateKey<DCRTPoly>;   // secret key of par.
using Public_Key = PublicKey<DCRTPoly>;     // public key of par.
using vecInt = std::vector<int64_t>;        // vector of intts

int main(int argc, char const *argv[]) {
	uint32_t multDepth = 2;
	uint32_t scaleFactorBits = 56;
	//uint32_t batchSize = 64; //power of two numbers: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512
	int i, j;
	int id = std::stoi(argv[1]);

	std::string fname;

	// Clinical dataset size
	int n_patients = 3691;
	int n_variables = 480; 
	
	SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
	CCParams<CryptoContextCKKSRNS> parameters;
    	parameters.SetMultiplicativeDepth(multDepth);
    	parameters.SetScalingModSize(scaleFactorBits);
	parameters.SetKeySwitchTechnique(HYBRID); 		//HYBRID
	parameters.SetScalingTechnique(FLEXIBLEAUTO); 	//FIXEDMANUAL FLEXIBLEAUTO
    	//parameters.SetBatchSize(batchSize);
	parameters.SetSecurityLevel(HEStd_128_classic);

	CC cc = GenCryptoContext(parameters);
	cc->Enable(PKE);
	cc->Enable(KEYSWITCH);
	cc->Enable(LEVELEDSHE); 	//EvalMultKeyGen
	cc->Enable(ADVANCEDSHE); 	//EvalSumKeyGen

	Key_Pair kp = cc->KeyGen();
	cc->EvalMultKeyGen(kp.secretKey);
	cc->EvalSumKeyGen(kp.secretKey);

	std::cout << "Batch size = " << cc->GetCryptoParameters()->GetEncodingParams()->GetBatchSize() << std::endl;

	std::vector<double> vector;
	PT pt;
	CT ct;
	std::vector<CT> ct_;

	std::random_device seed;
  	std::mt19937 gen{seed()}; // seed the generator
  	std::uniform_real_distribution<> dist{10,30};

	auto start = high_resolution_clock::now();
	for (i=0; i<n_patients; i++ ) {
        	vector.push_back(dist(gen));
        }
        pt = cc->MakeCKKSPackedPlaintext(vector);
        ct = cc->Encrypt(kp.publicKey, pt);
        ct_.push_back(ct);
	auto stop = high_resolution_clock::now();

	for (i=1; i<n_variables; i++ ) { 
		vector.clear();
		for (j=0; j<n_patients; j++ ) {
			vector.push_back(dist(gen));
		}
		pt = cc->MakeCKKSPackedPlaintext(vector);
		ct = cc->Encrypt(kp.publicKey, pt);
		ct_.push_back(ct);
	}

	auto duration = duration_cast<milliseconds>( stop - start );
	std::cout << "Time required to create one vector: " << duration.count() << " [ms]" << std::endl;

	start = high_resolution_clock::now();
	std::cout << "Serializing first vector to file" << std::endl;
    	if(!Serial::SerializeToFile("vector_0", ct_[0], SerType::BINARY)) {
      		std::cout << "Error writing to vector_0" << std::endl;
    	}
	stop = high_resolution_clock::now();
	auto duration2 = duration_cast<microseconds>( stop - start );
	std::cout << "Time required to serialize and store one vector: " << duration2.count() << " [us]" << std::endl;

	std::cout << "Serializing rest of the vectors" << std::endl;

        #pragma omp parallel for
        for (i=1; i<n_variables; i++) {
                fname = "vector_" + std::to_string(i);
                if(!Serial::SerializeToFile(fname, ct_[i], SerType::BINARY)) {
                         std::cout << "Error writing to " << fname << std::endl;
                }
        }
	
	std::vector<double> request, verification;

	// opt-out algorithm
	start = high_resolution_clock::now();
	for (i=0; i<n_patients; i++ ) {
		if (i == id) {
			request.push_back(0);
		}
		else {
			request.push_back(1);
		}
	}
	pt = cc->MakeCKKSPackedPlaintext(request);
	ct = cc->Encrypt(kp.publicKey, pt);
	
	std::vector<CT> opt(ct_.size());

	#pragma omp parallel for
	for (i=0; i<int(ct_.size()); i++) {
		//opt.push_back(cc->EvalMult(ct,ct_[i]));
		opt[i] = cc->EvalMult(ct,ct_[i]);
	}
	stop = high_resolution_clock::now();
	auto duration_optout = duration_cast<milliseconds>( stop - start );

	// Verification
	start = high_resolution_clock::now();
	for (i=0; i<n_patients; i++ ) {
                if (i == id) {
                        verification.push_back(1);
                }
                else {
                        verification.push_back(0);
                }
        }
	pt = cc->MakeCKKSPackedPlaintext(verification);
	ct = cc->Encrypt(kp.publicKey, pt);

	std::vector<CT> ver(opt.size());

	#pragma omp parallel for
	for (i=0; i<int(opt.size()); i++) {
                //ver.push_back(cc->EvalMult(ct,opt[i]));
                ver[i] = cc->EvalMult(ct,opt[i]);
        }
	stop = high_resolution_clock::now();
	auto duration_verification = duration_cast<milliseconds>( stop - start );

	// // Decryption
	Plaintext pt_proof;
	std::vector< std::vector<double> > proof;
	
	for (i=0; i<int(ver.size()); i++) {
		cc->Decrypt(kp.secretKey, ver[i], &pt_proof);
		proof.push_back(pt_proof->GetRealPackedValue());
	}

	// // Displaying the results
	std::vector<double> result;
	for (i=0; i<int(proof.size()); i++ ) {
		result.push_back(proof[i][id]);
	}
	std::cout << "Proof of deletion from the dataset: " << std::endl << result << std::endl << std::endl;
	std::cout << "Elapsed time: " << std::endl;
	std::cout << "Opt-out algorithm:     " << duration_optout.count() 	<< " [ms] " << std::endl;
	std::cout << "Verification algorithm " << duration_verification.count() << " [ms] " << std::endl;

	return 0;
}

