#include<cstdio>
#include<stdlib.h>
#include<omp.h>
#include<fstream>
#include<iostream>
#include<string>
#include<vector>
#include <unordered_map>

/*
    Implements the psuedolikelihood direct coupling analysis (plmDCA)
    for protein and RNA sequence families.

    Authors:
*/

void printSeqsBinaryState(const std::vector<std::vector<std::vector<unsigned int>>> & seqs_binary_form)
{
    /* 
    Prints elements of sequences that are in binary form. 
    This function is for testing purpose only. It will be romoved.
    */
    auto num_seqs = seqs_binary_form.size();
    auto seqs_len = seqs_binary_form[0].size();
    auto num_site_states = seqs_binary_form[0][0].size();
    std::cout << "***********num_seqs = " << num_seqs << std::endl; 
    std::cout << "***********num_site_states = " << num_site_states << std::endl;
    std::cout << "***********seqs_len = " << seqs_len << std::endl;
    std::cout << "******seqs_binary_form[0].size() = " << seqs_binary_form[0].size() << std::endl;
    std::cout << "******seqs_binary_form[0][0].size()=" << seqs_binary_form[0][0].size() << std::endl;
    for(unsigned int i = 0; i < num_seqs; ++i){
        std::cout << "sequence " << i + 1 << std::endl;
        for(unsigned int j = 0; j < seqs_len; ++j){
            std::cout << "site " << j + 1 << std::endl;
            for(unsigned int k = 0; k < num_site_states; ++k){
                std::cout << seqs_binary_form[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
    }
}

// obtain sequences length from FASTA file
unsigned int getSequencesLength(std::string msa_file)
{
    unsigned int seqs_len=0;
    std::ifstream msa_file_stream(msa_file);
    std::string current_line;
    int line_counter = 0;
    while(std::getline(msa_file_stream, current_line)){
        if(!current_line.empty()){
            if(current_line[0] != '>') ++line_counter;
            if(line_counter==1){ 
                seqs_len = current_line.length();
                break;
            }
        }
    }
    return seqs_len;
}


// put sequences in integer representation 
std::vector<std::vector<unsigned int>> getSequencesIntForm(const unsigned int biomolecule, const std::string  msa_file)
{
    std::unordered_map<char, unsigned int> res_mapping;
    if(biomolecule==1){
    /*  Protein residues mapping
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W':19, 'Y':20,
        '-':21, '.':21, '~':21,
    */
        res_mapping['A'] = 1;
        res_mapping['C'] = 2;
        res_mapping['D'] = 3;
        res_mapping['E'] = 4;
        res_mapping['F'] = 5;
        res_mapping['G'] = 6;
        res_mapping['H'] = 7;
        res_mapping['I'] = 8;
        res_mapping['K'] = 9;
        res_mapping['L'] = 10;
        res_mapping['M'] = 11;
        res_mapping['N'] = 12;
        res_mapping['P'] = 13;
        res_mapping['Q'] = 14;
        res_mapping['R'] = 15;
        res_mapping['S'] = 16;
        res_mapping['T'] = 17;
        res_mapping['V'] = 18;
        res_mapping['W'] = 19;
        res_mapping['Y'] = 20;
        res_mapping['-'] = 21;
        res_mapping['.'] = 21;
        res_mapping['~'] = 21;
    }else{
    /* RNA residues mapping
        'A':1, 'C':2, 'G':3, 'U':4, '-':5, '.':5, '~':5
    */
       res_mapping['A'] = 1;
       res_mapping['C'] = 2;
       res_mapping['G'] = 3;
       res_mapping['U'] = 4;
       res_mapping['-'] = 5;
       res_mapping['~'] = 5;
       res_mapping['.'] = 5;
    }
    std::vector<std::vector<unsigned int>> seqs_int_form;
    auto seqs_len = getSequencesLength(msa_file);
    std::ifstream msa_file_stream(msa_file);
    std::string current_line;
    std::string current_seq;
    int line_counter=0;
    std::vector<unsigned int> current_seq_int;

    while(std::getline(msa_file_stream, current_line)){
        if(!current_line.empty() && current_line[0] != '>'){ 
            for(unsigned int i = 0; i < seqs_len; ++i) {
                
                current_seq_int.emplace_back(res_mapping[ toupper(current_line[i])]);
            }
            // append the current sequences to seqs_int_form vector
            seqs_int_form.emplace_back(current_seq_int);
            // clear current_seq_int so that it is used to capture the next sequence
            current_seq_int.clear();
            // count number of sequences read from msa_file
            line_counter++;           
        }
    }
    std::cout << "Total number of sequences found: " << line_counter << std::endl;
    return seqs_int_form;
}


// Write sequences in binary form
std::vector<std::vector<std::vector<unsigned int>>> sequencesBinaryForm(const std::vector<std::vector<unsigned int>> & seqs_int_form, 
    const unsigned int num_site_states)
{
    std::vector<std::vector<std::vector<unsigned int>>> seqs_binary_form;
    auto seqs_len = seqs_int_form[0].size();
    auto num_seqs = seqs_int_form.size();
    std::vector< std::vector<unsigned int>> current_seq_binary_form;
    std::vector<unsigned int>  site_full_state(num_site_states);
    std::vector<unsigned int> current_seq;
    for (unsigned int i = 0; i < num_seqs; ++i){
        current_seq = seqs_int_form[i];
        for(unsigned int j = 0; j < seqs_len; ++j){
            for (unsigned int k = 0; k < num_site_states; ++k){
                site_full_state[k] = current_seq[j] == k + 1 ? 1 : 0;
            }
            current_seq_binary_form.emplace_back(site_full_state);
        }
        seqs_binary_form.emplace_back(current_seq_binary_form);
        current_seq_binary_form.clear();
    }
    seqs_binary_form.shrink_to_fit();
    return seqs_binary_form;
}

//compute sequences weigth
std::vector<double> computeSequencesWeight(const unsigned int biomolecule, std::string msa_file, double seqid)
{   
    /* 
    Computes the "normalized" weight of sequnces. 
    */
    auto seqs_int_form = getSequencesIntForm(biomolecule, msa_file);
    auto num_sequences = seqs_int_form.size();
    auto seqs_len = seqs_int_form[0].size();
    std::vector<double> seqs_weight(num_sequences);
    unsigned int num_identical_residues;
    
    //Assign all sequences of count = 1 
    for(unsigned int i=0; i < num_sequences; ++i){
        seqs_weight[i] = 1.0;
    }
    std::vector<unsigned int> seq_i, seq_j;
    double similarity_ij;
    for(unsigned int i=0; i < num_sequences - 1; ++i){
        seq_i = seqs_int_form[i];
        for(unsigned int j = i + 1; j < num_sequences; ++j){
            seq_j = seqs_int_form[j];
            num_identical_residues = 0;
            for(unsigned int site =0; site < seqs_len; ++site){
                if(seq_i[site] == seq_j[site]) num_identical_residues++;
            }
            similarity_ij = (double)num_identical_residues/(double)seqs_len;
            if (similarity_ij > seqid){
                seqs_weight[i] += 1.0;
                seqs_weight[j] += 1.0;
            }
        }
    }

    // "Normalize" sequences weight 
    for(unsigned int i=0; i < num_sequences; ++i){
        seqs_weight[i] = 1.0/seqs_weight[i];
    }
    return seqs_weight;
}


//compute effective number of sequences

double computeEffectiveNumSeqs(std::vector<double> const& seqs_weight)
{
    double eff_num_seqs = 0.0;
    auto num_seqs = seqs_weight.size();
    for(unsigned int i = 0 ; i< num_seqs; ++i) eff_num_seqs += seqs_weight[i];
    return eff_num_seqs;
}

//compute single site frequencies
std::vector<std::vector<double>> computeWeightedSingleSiteFreqs(std::vector<std::vector<unsigned int>> const& seqs_int_form, 
    std::vector<double> const& weights, const unsigned int num_site_states)
{
    std::vector<std::vector<double>> single_site_freqs;
    auto num_seqs = seqs_int_form.size();
    auto seqs_len = seqs_int_form[0].size();
    auto eff_num_seqs = computeEffectiveNumSeqs(weights);
    std::vector<double> current_site_freqs(num_site_states);
    double freq_ia;
    for(unsigned int i = 0; i < seqs_len; ++i){
        for(unsigned int j = 0; j < num_site_states; ++j){
            freq_ia = 0.0;
            for(unsigned int k = 0; k < num_seqs; ++k){
                if(seqs_int_form[k][i] == j + 1) freq_ia += weights[k];
            }
            freq_ia /= eff_num_seqs;
            current_site_freqs[j] = freq_ia;
        }
        single_site_freqs.emplace_back(current_site_freqs);
    }
    return single_site_freqs;
}

// Fields initializer
std::vector<std::vector<double>> initFields(const unsigned int seqs_len, const unsigned int num_site_states)
{
    std::vector<std::vector<double>> fields;
    std::vector<double> current_site_fields(num_site_states);
    for(unsigned int i = 0; i < num_site_states; ++i) current_site_fields[i] = 1.0/(double)num_site_states;
    for(unsigned int i = 0; i < seqs_len; ++i) fields.emplace_back(current_site_fields);
    fields.shrink_to_fit();
    return fields;

}

void printFields(const std::vector<std::vector<double>> & fields)
{
    auto seqs_len = fields.size();
    auto num_site_states = fields[0].size();
    for(unsigned int i = 0; i < seqs_len; ++i){
        std::cout << "site " << i + 1 << std::endl;
        for(unsigned int j = 0; j < num_site_states; ++j){
            std::cout << fields[i][j] << "," ;
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<std::vector<std::vector<double>>>> fourDimVecFactory(unsigned int vec_size_1, unsigned int vec_size_2, 
    unsigned int vec_size_3, unsigned int vec_size_4)
{
    std::vector<double> fourth_dim_vec(vec_size_4);
    std::vector<std::vector<double>> third_dim_vec(vec_size_3, fourth_dim_vec);
    std::vector<std::vector<std::vector<double>>> second_dim_vec(vec_size_2, third_dim_vec);
    std::vector<std::vector<std::vector<std::vector<double>>>> the_vector(vec_size_1, second_dim_vec);
    return the_vector;
}

// Couplings initializer 
std::vector<std::vector<std::vector<std::vector<double>>>> initCouplings(const unsigned int seqs_len, 
    const unsigned int num_site_states)
{  
   auto couplings = fourDimVecFactory(seqs_len, seqs_len, num_site_states, num_site_states);
    for(unsigned int i = 0; i < seqs_len; ++i){
        for(unsigned int j = 0 ; j < seqs_len; ++j){
            for(unsigned int a = 0; a < num_site_states; ++a){
                for(unsigned int b = 0; b < num_site_states; ++b){
                    couplings[i][j][a][b] =  i==j? 0.0 : (double)(num_site_states * num_site_states)/(double)seqs_len;
                }
            }
        }

    }
    return couplings; 
}


//print couplings
void printCouplings(std::vector<std::vector<std::vector<std::vector<double>>>> const& couplings)
{
    auto seqs_len = couplings[0].size();
    auto num_site_states = couplings[0][0][0].size();
    for(unsigned int i = 0; i < seqs_len - 1; ++i){
        for(unsigned int j = 0; j < seqs_len; ++j){
            for(unsigned int a = 0; a < num_site_states; ++a)
            {
                for(unsigned int b = 0;  b < num_site_states; ++b){
                    std::cout << i << " " << j << " " << a << " " << b << " " <<  couplings[i][j][a][b] << std::endl;
                }
            }
        }
    }
}

//compute weighted pair site frequencies 
std::vector<std::vector<std::vector<std::vector<double>>>> computeWeightedPairSiteFreqs(
    std::vector<std::vector<unsigned int>> const& seqs_int_form, 
    std::vector<double> const& weights, unsigned int num_site_states)
{
    auto eff_num_seqs = computeEffectiveNumSeqs(weights);
    auto num_seqs = seqs_int_form.size();
    auto seqs_len = seqs_int_form[0].size();
    auto pair_site_freqs = fourDimVecFactory(seqs_len, seqs_len, num_site_states, num_site_states);
    double freq_ij_ab;
    for(unsigned int i = 0; i < seqs_len; ++i){
        for(unsigned int j = 0; j < seqs_len; ++j){
            for(unsigned int a = 0;  a < num_site_states; ++a){
                for(unsigned int b = 0; b < num_site_states; ++b){
                    freq_ij_ab = 0.0;
                    for(unsigned int k = 0; k < num_seqs; ++k){
                        if(seqs_int_form[k][i] == a + 1 && seqs_int_form[k][j] == b + 1) freq_ij_ab += weights[k];
                    }
                    pair_site_freqs[i][j][a][b] = freq_ij_ab/eff_num_seqs;
                }
            }
        }
    }
    return pair_site_freqs;
}


//print pair-site frequencies
void printPairSiteFreqs(std::vector<std::vector<std::vector<std::vector<double>>>> const& pair_site_freqs)
{
    auto seqs_len = pair_site_freqs.size();
    auto num_site_states = pair_site_freqs[0][0][0].size();
}


void computeParams(std::vector<std::vector<double>>  & fields, 
    std::vector<std::vector<std::vector<std::vector<double>>>> & couplings)
{ 
    std::cout << "computing fileds and couplings" << std::endl;
}

extern "C" double* plmdcaBackend(const unsigned int biomolecule, unsigned int num_site_states, const char* msa_file, 
    const unsigned int seqs_len, const double seqid, const double lambda_h, const double lambda_J, int num_iter_steps)
{  
    /*
    The plmDCA main function. This function provides interface beteween the c++ 
    backend and the Python PlmDCA class. The arguments for this function are 
    attributes of the PlmDCA class.
    */ 
    std::cout << "***********--PLMDCA BACKEND--************" << std::endl;
    auto seqs_weight = computeSequencesWeight(biomolecule, msa_file, seqid);
    double eff_num_seqs = 0.0;
    for(unsigned int i=0; i < seqs_weight.size(); ++i) eff_num_seqs += seqs_weight[i];
    //for(unsigned int i = 0; i < seqs_weight.size(); ++i) std::cout << i << "\t" << seqs_weight[i] << std::endl;
    std::cout << "Effective number of sequences: " << eff_num_seqs << std::endl;
    //std::cout << "seqs_int_form vector size: " << seqs_int_form.size() << std::endl;
    std::cout << "plmDCA parameters" << std::endl;
    std::cout << "Sequences length:" << seqs_len << std::endl;
    std::cout << "lambda_h: " << lambda_h << std::endl;
    std::cout << "lamdba_J: " << lambda_J << std::endl;
    std::cout << "Number of site states: " << num_site_states << std::endl;
    std::cout << "Biomolecule: " << biomolecule << std::endl;
    auto seqs_int_form = getSequencesIntForm(biomolecule, msa_file);
    //auto seqs_binary_form = sequencesBinaryForm(seqs_int_form, num_site_states);
    //printSeqsBinaryState(seqs_binary_form);
    //auto fields = initFields(seqs_len, num_site_states);
    //std::cout << "Fields shape: (" << fields.size() << ", " << fields[0].size()<< ")" << std::endl;
    //printFields(fields);
    //computeWeightedSingleSiteFreqs(seqs_int_form, seqs_weight, num_site_states);
    //auto couplings = initCouplings(seqs_len, num_site_states);
    //printCouplings(couplings);
    computeWeightedPairSiteFreqs(seqs_int_form, seqs_weight, num_site_states);
    double * data_to_python = new double[seqs_len];
    for(unsigned int  i = 0; i < seqs_len; ++i) data_to_python[i] = 5.0 * i;
    std::cout << "***********--PLMDCA BACKEND--**************" << std::endl;
    return data_to_python;
}


int main()
{
    return 0;
}