
#if defined(_OPENMP)
    #include <omp.h>
#endif
#include <stdexcept> 
#include "include/plmdca.h"
#include <cstdio>

/*
    Implements pseudolikelihood maximization direct couplings analysis for 
    Protein and RNA sequence families.

    Author: Mehari B. Zerihun

*/

PlmDCA::PlmDCA(
    const char* m_msa_file, unsigned int m_biomolecule, unsigned int m_seqs_len, 
    unsigned int m_num_site_states, float m_seqid, float m_lambda_h, 
    float m_lambda_J, unsigned int m_num_threads
):msa_file(m_msa_file), biomolecule(m_biomolecule), seqs_len(m_seqs_len), 
    num_site_states(m_num_site_states), seqid(m_seqid), lambda_h(m_lambda_h),   
    lambda_J(m_lambda_J), num_threads(m_num_threads)
{
    /*  PlmDCA constructor

        Parameters
        ----------
            m_msa_file          : Path to FASTA formatted multiple sequence alignment file.
            m_biomolecule       : Type of biomolecule the MSA file represents.
            seqs_len            : Length of sequences in MSA. 
            m_num_site_states   : Number of possible states in a sequence of MSA.
            m_seqid             : Sequence identity threshold.
            m_lambda_h          : Value of fields regularization penality.
            m_lambda_J          : Value of couplings regularization penality.
            m_num_threads       : Number of threads (when OpenMP is supported)

    */

    this->num_fields = this->seqs_len * this->num_site_states;
    this->num_couplings = (this->seqs_len * (this->seqs_len -1))/2*(this->num_site_states * this->num_site_states);
    this->num_fields_and_couplings = this->num_fields + this->num_couplings;
    this->seqs_int_form = readSequencesFromFile();
    this->num_seqs = this->seqs_int_form.size();
    this->seqs_weight = this->computeSeqsWeight();
    this->Meff = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.f);
    
}


std::vector<std::vector<float>> PlmDCA::getSingleSiteFreqs()
{
    /*Computes single site frequencies.

    Parameters
    ----------
        this    : An instance of PlmDCA class
    
    Returns
    -------
        single_site_reqs    : Normalized single site frequencies.
    */
    std::vector<std::vector<float>> single_site_freqs(this->seqs_len);
    std::vector<float> current_site_freqs(this->num_site_states);
    for(unsigned int i = 0; i < this->seqs_len; ++i){
        std::fill(current_site_freqs.begin(), current_site_freqs.end(), 0.f);
        for(unsigned int n = 0; n < this->num_seqs; ++n){
            auto a = this->seqs_int_form[n][i];
            current_site_freqs[a] += this->seqs_weight[n]; 
        }
        single_site_freqs[i] = current_site_freqs;
    }
    // Divide frequency counts by effective number of sequences

    for(unsigned int i = 0; i < this->seqs_len; ++i){
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            single_site_freqs[i][a] /= this->Meff;
        }
    }
    return single_site_freqs; 
}


//Compute pair-site frequencies

std::vector<float> PlmDCA::getPairSiteFreqs()
{
    /*Computes pair-site frequencies

    Parameters
    ----------
        this    : PlmDCA instance

    Returns
    -------
        fij     : Pair-site frequencies

    */
    auto const& L = this->seqs_len;
    auto const& q = this->num_site_states;
    auto const& Nseq = this->num_seqs;
    auto const freqs_size = (L*(L - 1)/2)*q*q;
    std::vector<float> fij(freqs_size, 0.f);
    float Meff = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.f);
    float MeffInv = 1.0/Meff;
    
    #pragma omp parallel for num_threads(this->num_threads)
    for(unsigned int i = 0; i < L; ++i){
        std::vector<float> current_fij(L*q*q, 0.f);
        for(unsigned int n = 0; n < Nseq; ++n){
            auto const& current_seq_weight = this->seqs_weight[n];
            auto const& current_seq = this->seqs_int_form[n];
            for(unsigned int j = 0;  j < i; ++j){
                current_fij[this->mapIndexPairSiteFreqsLocal(j, current_seq[j], current_seq[i])] += current_seq_weight;
            }
            for(unsigned int j = i + 1; j <  L; ++j){
                current_fij[this->mapIndexPairSiteFreqsLocal(j, current_seq[i], current_seq[j])] += current_seq_weight;
            }
        }//iteration over sequences

        #pragma omp critical
        {
        for(unsigned int j = 0; j < i; ++j){
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    fij[this->mapIndexPairSiteFreqs(j, i, a, b)] =  MeffInv * current_fij[this->mapIndexPairSiteFreqsLocal(j, a, b)];
                }
            }
        }
        for(unsigned int j = i + 1; j < L; ++j){
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    fij[this->mapIndexPairSiteFreqs(i, j, a, b)] =  MeffInv * current_fij[this->mapIndexPairSiteFreqsLocal(j, a, b)];
                }
            }
        }
        }//end of critical section
    }// iteration over sites
    return fij; 
}


//Compute pair-site frequencies and store them in a 4d vector
std::vector<std::vector<std::vector<std::vector<float>>>> PlmDCA::getPairSiteFreqsFragmented()
{
    /*Compute pair-site frequencies and store them in a 4d vector.
    This is used for comparing with the multithreaded implementation version,
    otherwise it will not be used in production.
    */

    float eff_num_seqs = std::accumulate(this->seqs_weight.begin(),  this->seqs_weight.end(), 0.f);
    auto const& weights = this->seqs_weight;
    auto const& num_seqs = this->num_seqs;
    auto const& seqs_len = this->seqs_len;
    auto pair_site_freqs = this->fourDimVecFactory(this->seqs_len, this->seqs_len, 
        this->num_site_states, this->num_site_states
    );
    float freq_ij_ab;
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


// 4d vector factory 
std::vector<std::vector<std::vector<std::vector<float>>>> PlmDCA:: fourDimVecFactory(
    unsigned int const vec_size_1, unsigned int const vec_size_2, 
    unsigned int const vec_size_3, unsigned int const vec_size_4
)
{
    std::vector<float> fourth_dim_vec(vec_size_4);
    std::vector<std::vector<float>> third_dim_vec(vec_size_3, fourth_dim_vec);
    std::vector<std::vector<std::vector<float>>> second_dim_vec(vec_size_2, third_dim_vec);
    std::vector<std::vector<std::vector<std::vector<float>>>> the_vector(vec_size_1, second_dim_vec);
    return the_vector;
}


// Investigate single site frequencies
void PlmDCA::testSingleSiteFreqs()
{
    auto ssfreqs = this->getSingleSiteFreqs();
    float Meff = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.f);
    std::cout << "Effective number of sequences = " << Meff << std::endl;
    for(unsigned int i = 0;  i < this->seqs_len; ++i){
        float local_freqs_sum = 0.0;
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            local_freqs_sum += ssfreqs[i][a];
        }
        //std::cout << i << " sum of freqs = " << local_freqs_sum << std::endl;
    }
}


// Initialize fields and couplings
void PlmDCA::initFieldsAndCouplings(float* fields_and_couplings)
{
    /*Sets initial value for fields and couplings

    Parameters
    ----------
        this                    : An instance of PlmDCA class.
        fields_and_couplings    : Array of fields and couplings.

    Returns
    -------
        void    : No return value

    */
        auto single_site_freqs =  this->getSingleSiteFreqs();
        unsigned int index;
        float epsilon = 1.f;
        for(unsigned int i = 0; i < this->seqs_len; ++i){
            for(unsigned int a = 0; a < this->num_site_states; ++a){
                index = a + this->num_site_states * i;
                fields_and_couplings[index] = std::log(single_site_freqs[i][a] * this->Meff  + epsilon);
            }
        }

        for(unsigned int i = 0; i < this->seqs_len; ++i){

            auto start_indx = i * this->num_site_states;
            auto end_indx = start_indx + this->num_site_states;
            auto hi_sum = std::accumulate(fields_and_couplings + start_indx, fields_and_couplings + end_indx, 0.f);
            auto hi_sum_av = hi_sum/this->num_site_states;
            for(unsigned int a = 0; a < this->num_site_states; ++a){
                auto index = a + this->num_site_states * i;
                fields_and_couplings[index] -= hi_sum_av;
            }
        }

        
        for(unsigned int i = this->num_fields; i < this->num_fields_and_couplings; ++i){
            fields_and_couplings[i] = 0.f;

        }
    
}


//Map pair-site frequencies indices
unsigned int PlmDCA::mapIndexPairSiteFreqs(const unsigned int i, const unsigned int j, 
    const unsigned int a, const unsigned int b)
{
    /*Maps pair site indices into contiguous memory layout. The mapping is an upper
    triangular matrix form with i < j. Each site pair (i, j) occupies the q*q float 
    consecutive memory block.

        Parameters
        ----------
            this    : PlmDCA instance
            i       : Site i in site pair (i, j) with j > i.
            j       : Site j in site pair (i, j) with j > i.
            a       : Residue at site i. 
            b       : Residue at site j.

        Returns
        -------
            k   : Site location in contiguous memory layout for site-pair (i, j)
            and corresponding residue pair (a, b)  
    */
    auto q = this->num_site_states;
    auto L = this->seqs_len; 
    auto site = ((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q;
    auto k = site + b + a * q;
    return k;
}


//Map indices of pair-site frequencies for one site paired with all others
unsigned int PlmDCA::mapIndexPairSiteFreqsLocal(const unsigned int j, 
    const unsigned int a, const unsigned int b)
{
    /*Maps local pair-site frequencies indices

    Parameters
    ----------
        this    : An instance of PlmDCA class.
        j       : Outer most loop index.
        a       : Residue at site i when for site (i, j) with j > i.
        b       : Residue at site j.
    */
   auto const&  q = this->num_site_states;
   return  b +  q * ( a + q * j);

}


void PlmDCA::printMapIndexPairSiteFreqsLocal(const unsigned int i)
{
    for(unsigned int j = 0; j <  i; ++j){
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            for(unsigned int b = 0;  b < this->num_site_states;++b){
                std::cout << j << " " << a << " " << b << " " << this->mapIndexPairSiteFreqsLocal(j , a, b) << std::endl;
            }
        }
    }
    for(unsigned int j = i + 1; j <  this->seqs_len; ++j){
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            for(unsigned int b = 0;  b < this->num_site_states; ++b){
                std::cout << j << " " << a << " " << b << " " << this->mapIndexPairSiteFreqsLocal(j , a, b) << std::endl;
            }
        }
    }
}

// Map couplings indices from upper triangular matrix to one-dimensional vector
unsigned int PlmDCA::mapIndexCouplings(const unsigned int i, const unsigned int j, 
    const unsigned int a, const unsigned int b)
{   
    /*  Couplings index mapper.

        Parameters
        ----------
            this    : An instance of PlmDCA class.
            i       : Site i for site pair (i, j) such that j > i. 
            j       : Site j for site pair (i, j) such that j > i.
            a       : Residue/gap at the trailing site i
            b       : Residue/gap at site j
        
        Returns
        -------
            k   : Index mapping from upper triangular matrix representation to 
                that of one-dimensional vector
    */

    auto q = this->num_site_states;
    auto L = this->seqs_len; 
    auto site = ((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q;
    auto  k =  L * q + site + b + a * q;
    return k;
}


// Map fields index from a two dimensional to one-dimensional vector
unsigned int PlmDCA::mapIndexFields(const unsigned int i, const unsigned int a)
{
    /*  Fields index mapper.

        Parameters
        ----------
            this    : An instance of PlmDCA class. 
            i       : Site in sequence of MSA
            a       : Residue/gap at site i

        Returns
        --------
            k   : Index mapping from two-dimensional representation of fileds to 
                one dimensional representation.
    */
    auto q =  this->num_site_states;
    unsigned int k = a + i * q; 
    return k; 
}


// Map couplings index per site
unsigned int PlmDCA::mapIndexCouplingsOneSite(const unsigned int j, const unsigned int a, 
    const unsigned int b)
{
    /* Mapps couplings when a particular site is coupled with all other sites in
    sequence. 

    Parameters
    ----------
        j   : Looping index that pairs a particular site with all other remaining
            sites.
        a   : Residue at site i for site pair (i, j) such that j > i.
        b   : Residue at site j.

    Returns
    -------
        k : Mapped site index 
    */
    auto const&  q = this->num_site_states;
    unsigned int k = b + q * ( a +  q * j);
    return k;
}


void PlmDCA::printSeqs()
{
    auto M = this->seqs_int_form.size();
    auto L = this->seqs_int_form[0].size();
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    
    std::cout << "#num seqs int form:    " << M << std::endl;
    std::cout << "#len seqs:             " << L << std::endl;
    for(unsigned int i = 0; i < M; ++i){
        for(unsigned int j = 0; j < L; ++j){
            //std::cout << seqs_int_form[i][j] << ',';
        }
       // std::cout << std::endl;
    }
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
}


void PlmDCA::printIndexMappingFields()
{
    for(unsigned int i  = 0; i < this->seqs_len; ++i){
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            std::cout << i << " " << a << " " << this->mapIndexFields(i, a) << std::endl;
        }
    }
}


void PlmDCA::printIndexMappingCouplings()
{
    for(unsigned int i = 0; i <  this->seqs_len; ++i){
        for(unsigned int j =  i + 1; j < this->seqs_len; ++j){
            for(unsigned int a = 0; a < this->num_site_states; ++a){
                for(unsigned int b = 0; b < this->num_site_states; ++b){
                    std::cout <<  i << " " << " " << j <<  " " << a  << " " << 
                        b << " " << this->mapIndexCouplings(i, j, a, b) << std::endl;
                 }
            }
        }
    }
}


//compute gradients from site probabilities
float  PlmDCA::gradient(const float* fields_and_couplings, float* grad)
{
    /*Computes the gradient of the negative psuedolikelihood from alignment data 
    plus contributions from L2-norm regularization for the fields and couplings.

    Parameters
    ----------
        this                    : An instance of PlmDCA class
        fields_and_couplings    : Array of fields and couplings
        grad                    : Array of gradients 

    Returns
    -------
        fx                      : Value of objective function

    */
   

    auto const& q = this->num_site_states;
    auto const& L = this->seqs_len;
    auto const& Nseq = this->num_seqs;
    auto const& lh = this->lambda_h;
    auto const& lJ = this->lambda_J;

    float fx = 0.f;
    
    // Gradients of L2-norm regularization terms.
    for(unsigned int i = 0; i < L; ++i){
        for(unsigned int a = 0; a < q; ++a){
            //auto const& indx_ia = this->mapIndexFields(i, a);
            //auto const& hia = fields_and_couplings[indx_ia];
            auto const hia = fields_and_couplings[a + i * q];
            grad[a + q * i] = 2.f * lh * hia;
            fx += lh *  hia * hia;
        }
    }

    for(unsigned int i = 0; i < L - 1; ++i){
        for(unsigned int j = i + 1; j < L; ++j){
            auto k = L * q + ((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q;
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    //auto const& indx_ij_ab = this->mapIndexCouplings(i, j, a, b);
                    //auto const& Jijab = fields_and_couplings[indx_ij_ab];
                    auto const& Jijab = fields_and_couplings[k + b + q * a];
                    grad[k + b + q * a] = 2.f * lJ * Jijab;
                    fx += lJ *  Jijab * Jijab;
                }
            }
        }
    }


    // Compute gradients of the negative of pseudolikelihood from alignment data. 
    #pragma omp parallel for num_threads(this->num_threads)
    for(unsigned int i = 0; i < L; ++i){
        std::vector<float> prob_ni(q, 0.f);
        std::vector<float> fields_gradient(q, 0.f);
        std::vector<float> couplings_gradient(L * q * q, 0.f);
        float  fxi = 0.f;
        for(unsigned int n = 0; n < Nseq; ++n){
            // compute probability at site i for sequence n
            auto const& current_seq = this->seqs_int_form[n];
            for(unsigned int a = 0; a < q; ++a) prob_ni[a] += fields_and_couplings[a + i * q];            
            
            for(unsigned int j = 0; j < i; ++j){
                auto k = L * q + ((L *  (L - 1)/2) - (L - j) * ((L-j)-1)/2  + i  - j - 1) * q * q;
                for(unsigned int a = 0; a < q; ++a){
                    // Index mapping hint
                    // auto const& indx_jia = this->mapIndexCouplings(j, i, current_seq[j], a);
                    prob_ni[a] += fields_and_couplings[k + a + current_seq[j] * q];
                }
            }

            for(unsigned int j = i + 1; j < L; ++j){
                auto k = L * q + ((L * (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q;
                for(unsigned int a = 0; a < q; ++a){
                    // Index mapping hint
                    // auto const&  indx_ija = this->mapIndexCouplings(i, j, a, current_seq[j]);
                    prob_ni[a] += fields_and_couplings[k + current_seq[j]  + a * q];
                }
            }

            auto max_exp = prob_ni[0];
            for(unsigned int a = 0; a < q; ++a){
                if(prob_ni[a] > max_exp) {
                    max_exp = prob_ni[a];
                }
            } 
            for(unsigned int a = 0; a <  q; ++a) prob_ni[a] = std::exp(prob_ni[a] - max_exp);
            
            float zni = 0.f;
            for(unsigned int a = 0; a < q; ++a) zni += prob_ni[a];
            zni = 1.f/zni;
            for(unsigned int a = 0; a < q; ++a) prob_ni[a] *= zni;
            
            //compute the gradient of the minus log likelihood with respect to fields and couplings
            auto current_seq_weight = this->seqs_weight[n];
            auto const& res_i = current_seq[i];

            fxi -= current_seq_weight * std::log(prob_ni[res_i]);

            fields_gradient[res_i] -= current_seq_weight;
            for(unsigned int a = 0; a < q; ++a) fields_gradient[a] += current_seq_weight * prob_ni[a];
            
            for(unsigned int j = 0; j < i; ++j){
                //auto const& indx_ji = this->mapIndexCouplingsOneSite(j, current_seq[j], res_i);
                auto k = res_i + q * (current_seq[j] + q * j);
                couplings_gradient[k] -= current_seq_weight;
            }
            for(unsigned int j = i + 1; j < L; ++j){
                // auto const& indx_ij = this->mapIndexCouplingsOneSite(j, res_i,  current_seq[j]);
                auto k = current_seq[j] + q * (res_i +  q * j);
                couplings_gradient[k] -= current_seq_weight;
            }

            for(unsigned int j = 0; j < i; ++j){
                //auto k = b + q * ( a +  q * j);
                auto k =  q * (current_seq[j] +  q * j);
                for(unsigned int a = 0; a < q; ++a){
                    //auto const& indx_jia = this->mapIndexCouplingsOneSite(j, current_seq[j], a);
                    couplings_gradient[a + k] += current_seq_weight * prob_ni[a]; 
                }
            }
            for(unsigned int j = i + 1; j < L; ++j){
                //auto k = b + q * ( a +  q * j);
                for(unsigned int a = 0; a < q; ++a){
                    //auto const& indx_ija = this->mapIndexCouplingsOneSite(j, a, current_seq[j]);
                    auto k = current_seq[j] + q * (a +  q * j);
                    couplings_gradient[k] += current_seq_weight * prob_ni[a];
                }
            }   
        }// End of iteration through sequences
        
        #pragma omp critical
        {
            fx += fxi;
            for(unsigned int a = 0; a < q; ++a){
                //auto const& indx_ia = this->mapIndexFields(i, a);
                //grad[indx_ia] += fields_gradient[a];
                grad[a + q * i] += fields_gradient[a];

            }

            for(unsigned int j = 0; j < i; ++j){
                auto k = L * q + ((L *  (L - 1)/2) - (L - j) * ((L-j)-1)/2  + i  - j - 1) * q * q;
                for(unsigned int a = 0; a < q; ++a){
                    auto k_2 = q * (a + q * j);
                    for(unsigned int b = 0; b < q; ++b){
                        //auto const& indx_jiab = this->mapIndexCouplings(j, i, a, b);
                        //grad[indx_jiab] += couplings_gradient[this->mapIndexCouplingsOneSite(j, a, b)];
                        grad[k + b + a * q] += couplings_gradient[b + k_2];
                    }
                }
            }
            for(unsigned int j = i + 1; j < L; ++j){
                auto k = L * q + ((L *  (L - 1)/2) - (L - i) * ((L-i)-1)/2  + j  - i - 1) * q * q;
                for(unsigned int a = 0; a < q; ++a){
                    auto k_2 = q * (a + q * j);
                    for(unsigned int b = 0; b < q; ++b){
                        //auto const& indx_ijab = this->mapIndexCouplings(i, j, a, b);
                        //grad[indx_ijab] += couplings_gradient[this->mapIndexCouplingsOneSite(j, a, b)];
                        grad[k + b + a * q] += couplings_gradient[b + k_2];
                    }
                }
            }
        }//end of critical block
    
    } // End of iteration through sites

    return fx;
 }


//compute sequences weight
std::vector<float> PlmDCA::computeSeqsWeight()
{
    /*  Computes sequences weight.

        Parameters
        ----------
            this            : An instance of PlmDCA class.

        Returns
        -------
            m_seqs_weight   : The weigh of sequences computed with respect to 
                sequence identity obtained from this->seqid.

    */
    std::vector<float> m_seqs_weight(this->num_seqs);

    #if defined(_OPENMP)
        //initialize weights to zero for the parallel case since each sequences 
        //is going to be compared with itself.
        for(unsigned int i=0; i < this->num_seqs; ++i){
            m_seqs_weight[i] = 0.f; 
        }
        #pragma omp parallel for num_threads(this->num_threads)
        for(unsigned int i = 0; i < this->num_seqs; ++i){
            float similarity_ij;
            unsigned int num_identical_residues;
            for(unsigned int j = 0; j < this->num_seqs; ++j){
                num_identical_residues = 0;
                for(unsigned int site =0; site < this->seqs_len; ++site){
                    if(this->seqs_int_form[i][site] == this->seqs_int_form[j][site]) num_identical_residues++;
                }
                similarity_ij = (float)num_identical_residues/(float)this->seqs_len;
                if (similarity_ij > this->seqid) m_seqs_weight[i] += 1.f;
            }
        }
    #else
        // Initialize weights to one for unique sequence pair comparisons.
        for(unsigned int i = 0; i < this->num_seqs; ++i){
            m_seqs_weight[i] = 1.f;
        }
        float similarity_ij;
        unsigned int num_identical_residues;
        for(unsigned int i = 0; i < this->num_seqs - 1; ++i){
            for(unsigned int j = i + 1; j < this->num_seqs; ++j){
                num_identical_residues = 0;
                for(unsigned int site =0; site < this->seqs_len; ++site){
                    if(this->seqs_int_form[i][site] == this->seqs_int_form[j][site]) num_identical_residues++;
                }
                similarity_ij = (float)num_identical_residues/(float)this->seqs_len;
                if (similarity_ij > this->seqid){
                    m_seqs_weight[i] += 1.f;
                    m_seqs_weight[j] += 1.f;
                }
            }
        }
    #endif
    
    //"Normalize" sequences weight 
    for(unsigned int i = 0; i < this->num_seqs; ++i) m_seqs_weight[i] = 1.f/m_seqs_weight[i];
    return m_seqs_weight;
}



//print sequences for debugging 
void PlmDCA::printWeights()
{
    for(unsigned int i = 0; i < this->num_seqs; ++i){
        std::cout << "# " << i + 1 << " " << this->seqs_weight[i] << std::endl;
    }
    std::cout << "#Meff = " << std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.f) << std::endl;
}

//Read alignment from FASTA file
std::vector<std::vector<unsigned int>>  PlmDCA::readSequencesFromFile()
{
    /*  Reads sequences from FASTA file and mappes them into integer representation.
        All residues other than the standard residues are mapped to gap state value.

        Parameters
        ----------
            this    : An instance of PlmDCA class. 
        
        Returns
        -------
            seqs_int_form   : Integer represenation of sequences read from FASATA
                formatted alignment file. 
    */
    std::unordered_map<char, unsigned int> res_mapping;
    if(this->biomolecule==1){
    /*  Protein residues mapping (value minus 1 for optimization)
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W':19, 'Y':20,
        '-':21, '.':21, '~':21,
    */
        res_mapping['A'] = 0; res_mapping['C'] = 1; res_mapping['D'] = 2;
        res_mapping['E'] = 3; res_mapping['F'] = 4; res_mapping['G'] = 5;
        res_mapping['H'] = 6; res_mapping['I'] = 7; res_mapping['K'] = 8;
        res_mapping['L'] = 9; res_mapping['M'] = 10; res_mapping['N'] = 11;
        res_mapping['P'] = 12; res_mapping['Q'] = 13; res_mapping['R'] = 14;
        res_mapping['S'] = 15; res_mapping['T'] = 16; res_mapping['V'] = 17;
        res_mapping['W'] = 18; res_mapping['Y'] = 19; res_mapping['-'] = 20;
        res_mapping['.'] = 20; res_mapping['~'] = 20; res_mapping['B'] = 20;
        res_mapping['J'] = 20; res_mapping['O'] = 20; res_mapping['U'] = 20;
        res_mapping['X'] = 20; res_mapping['Z'] = 20;
    }else{
    /* RNA residues mapping
        'A':1, 'C':2, 'G':3, 'U':4, '-':5, '.':5, '~':5
    */
       res_mapping['A'] = 0; res_mapping['C'] = 1; res_mapping['G'] = 2;
       res_mapping['U'] = 3; res_mapping['-'] = 4; res_mapping['~'] = 4;
       res_mapping['.'] = 4; res_mapping['B'] = 4; res_mapping['D'] = 4;
       res_mapping['E'] = 4; res_mapping['F'] = 4; res_mapping['H'] = 4;
       res_mapping['I'] = 4; res_mapping['J'] = 4; res_mapping['K'] = 4;
       res_mapping['L'] = 4; res_mapping['M'] = 4; res_mapping['N'] = 4;
       res_mapping['O'] = 4; res_mapping['P'] = 4; res_mapping['Q'] = 4;
       res_mapping['R'] = 4; res_mapping['S'] = 4; res_mapping['T'] = 4;
       res_mapping['V'] = 4; res_mapping['W'] = 4; res_mapping['X'] = 4;
       res_mapping['Y'] = 4; res_mapping['Z'] = 4;
    }
    //TODO obtain the number of sequences and reserver memory to seqs_int_form
    //to avoid reallocation of memorry.
    std::vector<std::vector<unsigned int>> seqs_int_form;
    std::ifstream msa_file_stream(this->msa_file);
    std::string current_line;
    std::string current_seq;
    int line_counter=0;
    int unique_seq_counter = 0;
    std::vector<unsigned int> current_seq_int;

    if(msa_file_stream.fail()){
        std::cerr << "Unable to open file " << this->msa_file << std::endl;
        throw std::runtime_error("Unable to open file containing the MSA data\n");
    }

    while(std::getline(msa_file_stream, current_line)){
        if(!current_line.empty() && current_line[0] != '>'){ 
            for(unsigned int i = 0; i < seqs_len; ++i){
                
                current_seq_int.emplace_back(res_mapping.at(toupper(current_line[i])));
            }
    
            // take only unique sequences 
            if(std::find(seqs_int_form.begin(), seqs_int_form.end(), current_seq_int) == seqs_int_form.end()){
                seqs_int_form.emplace_back(current_seq_int);
                ++unique_seq_counter;
            }
            // clear current_seq_int so that it is used to capture the next sequence
            current_seq_int.clear();
            // count number of sequences read from msa_file
            line_counter++;           
        }
    }
    return seqs_int_form;
}
