#ifndef PLMDCA_BACKEND_H
#define PLMDCA_BACKEND_H
    #include<cstdio>    
    #include<cstdlib>
    #include<fstream>
    #include<iostream>
    #include<string>
    #include<vector>
    #include<unordered_map>
    #include<algorithm>
    #include<cmath>
    #include<numeric>
    #include"../lbfgs/include/lbfgs.h"


    class PlmDCA {

        public:
            
            
            PlmDCA(const char* m_msa_file, unsigned int m_biomolecule, unsigned int m_seqs_len, 
                unsigned int m_num_site_states, float m_seqid, float m_lambda_h, 
                float m_lambda_J, unsigned int num_threads
            );

            //float value(const TVector &fields_and_couplings);
            float gradient(const float* fields_and_couplings, float* grad);

            void initFieldsAndCouplings(float* fields_and_couplings);
            std::vector<std::vector<float>> getSingleSiteFreqs();
            std::vector<float> getPairSiteFreqs();
            std::vector<std::vector<std::vector<std::vector<float>>>> getPairSiteFreqsFragmented();
            std::vector<std::vector<std::vector<std::vector<float>>>> fourDimVecFactory(
                unsigned int const vec_size_1, unsigned int const vec_size_2, 
                unsigned int const vec_size_3, unsigned int const vec_size_4
            );
            void printPairSiteFreqsMultiThreaded();
            void printPairSiteFreqsFragmented();
            unsigned int mapIndexPairSiteFreqs(const unsigned int, const unsigned int,  
                const unsigned int, const unsigned int
            );
            unsigned int mapIndexPairSiteFreqsLocal(const unsigned int, const unsigned int,
                const unsigned int
            );
            void printMapIndexPairSiteFreqsLocal(const unsigned int);
            void testSingleSiteFreqs();
            std::vector<std::vector<unsigned int>> readSequencesFromFile();
            unsigned int mapIndexCouplings(const unsigned  int i, const unsigned int j, 
                const unsigned int a, const unsigned int b
            );
            unsigned int mapIndexCouplingsOneSite(const unsigned int j, 
                const unsigned int a, const unsigned int b
            );
            unsigned int mapIndexFields(const unsigned int i, const unsigned int a);
            void printIndexMappingFields();

            void printIndexMappingCouplings();

            std::vector<float> computeSeqsWeight();
            void printWeights();
            void runPlmDCALocal(unsigned int num_iteration);
            float* computeGradient();
            //void printDCAScores(float* h_and_J);
            int runPlmDCA();

            void printSeqs();
            
        private:
            const char* msa_file;
            unsigned int biomolecule;
            unsigned int seqs_len;
            unsigned int  num_site_states;
            float seqid;
            float lambda_h;
            float lambda_J;
            std::vector<std::vector<unsigned int>> seqs_int_form;
            unsigned int num_seqs;
            std::vector<float> seqs_weight;
            unsigned int  num_fields;
            unsigned int num_couplings;
            unsigned int num_fields_and_couplings;
            unsigned int  num_threads;
            float  Meff;
    };
#endif
