
#if defined(_OPENMP)
    #include <omp.h>
#endif
#include <stdexcept> 
#include "plmdca.h"
#include <stdio.h>

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

        Attributes
        ----------
            msa_file        : Path to FASTA formatted multiple sequence alignment file.
            biomolecule     : Type of biomolecule the MSA file represents.
            seqs_len        : Length of sequences in MSA. 
            num_site_states : Number of possible states in a sequence of MSA.
            seqid           : Sequence identity threshold.
            lambda_h        : Value of fields regularization penality.
            lambda_J        : Value of couplings regularization penality.    

    */
    this->num_fields = this->seqs_len * this->num_site_states;
    this->num_couplings = (this->seqs_len * (this->seqs_len -1))/2*(this->num_site_states * this->num_site_states);
    this->num_fields_and_couplings = this->num_fields + this->num_couplings;
    this->seqs_int_form = readSequencesFromFile();
    this->num_seqs = this->seqs_int_form.size();
    this->seqs_weight = this->computeSeqsWeight();
    
    /*
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "#msa file name:    " << this->msa_file << std::endl;
    std::cout << "#biomolecule:      " << this->biomolecule << std::endl;
    std::cout << "#sequences length: " << this->seqs_len << std::endl;
    std::cout << "#num site states:  " << this->num_site_states << std::endl;
    std::cout << "#sequence identity:" << this->seqid << std::endl;
    std::cout << "#num seqs:        " << this->num_seqs << std::endl;
    std::cout << "#lambda_h:        " << this->lambda_h << std::endl;
    std::cout << "#lamdba_J:        " << this->lambda_J << std::endl;
    std::cout << "#num_fields:      " << this->num_fields << std::endl;
    std::cout << "#num_couplings:   " << this->num_couplings << std::endl; 
    std::cout << "#num_fields_and_couplings:    " << this->num_fields_and_couplings << std::endl;
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    */
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
        std::fill(current_site_freqs.begin(), current_site_freqs.end(), 0.0);
        for(unsigned int n = 0; n < this->num_seqs; ++n){
            auto a = this->seqs_int_form[n][i] - 1;
            current_site_freqs[a] += this->seqs_weight[n]; 
        }
        single_site_freqs[i] = current_site_freqs;
    }
    // Divide frequency counts by effective number of sequences
    float eff_num_seqs = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.0);
    for(unsigned int i = 0; i < this->seqs_len; ++i){
        for(unsigned int a = 0; a < this->num_site_states; ++a){
            single_site_freqs[i][a] /= eff_num_seqs;
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
    std::vector<float> fij(freqs_size, 0.0);
    float Meff = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.0);
    float MeffInv = 1.0/Meff;
    
    #pragma omp parallel for num_threads(this->num_threads)
    for(unsigned int i = 0; i < L; ++i){
        std::vector<float> current_fij(L*q*q, 0.0);
        for(unsigned int n = 0; n < Nseq; ++n){
            auto const& current_seq_weight = this->seqs_weight[n];
            auto const& current_seq = this->seqs_int_form[n];
            for(unsigned int j = 0;  j < i; ++j){
                current_fij[this->mapIndexPairSiteFreqsLocal(j, current_seq[j] - 1, current_seq[i] - 1)] += current_seq_weight;
            }
            for(unsigned int j = i + 1; j <  L; ++j){
                current_fij[this->mapIndexPairSiteFreqsLocal(j, current_seq[i] - 1, current_seq[j] - 1)] += current_seq_weight;
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

    float eff_num_seqs = std::accumulate(this->seqs_weight.begin(),  this->seqs_weight.end(), 0.0);
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



//Print pair-site frequencies computed using the multithreaded version of pair-site
//frequency computation method
void PlmDCA::printPairSiteFreqsMultiThreaded()
{
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    auto fij = this->getPairSiteFreqs();
    auto const& L = this->seqs_len;
    auto const& q = this->num_site_states;
    for(unsigned int i = 0; i < L - 1; ++i){
        for(unsigned int j = i + 1; j < L; ++j){
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    auto indx = this->mapIndexPairSiteFreqs(i, j, a, b);
                    std::cout << i << " " << j << " " << a << " " << " " << b << " " << indx << " " << fij.at(indx) << std::endl;
                 }
            }
        }
    }
}


//Print pair-site frequencies computed using the fragmented memory layout method
void PlmDCA::printPairSiteFreqsFragmented()
{
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    auto fij_fragmented = this->getPairSiteFreqsFragmented();
    auto const& L = this->seqs_len;
    auto const& q = this->num_site_states;
    for(unsigned int i = 0; i < L - 1; ++i){
        for(unsigned int j = i + 1; j < L; ++j){
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    auto indx = this->mapIndexPairSiteFreqs(i, j, a, b);
                    std::cout << i << " " << j << " " << a << " " << " " << b << " " << indx << " " << fij_fragmented[i][j][a][b] << std::endl;
                 }
            }
        }
    }
}




// Investigate single site frequencies
void PlmDCA::testSingleSiteFreqs()
{
    auto ssfreqs = this->getSingleSiteFreqs();
    float Meff = std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.0);
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
        auto single_site_freqs =  this->getSingleSiteFreqs();
        unsigned int index;
        float epsilon = 0.00001;
        for(unsigned int i = 0; i < this->seqs_len; ++i){
            for(unsigned int a = 0; a < this->num_site_states; ++a){
                index = a + this->num_site_states*i;
                fields_and_couplings[index] = std::log(single_site_freqs[i][a] + epsilon);
                //std::cout << index << " " << this->fields_and_couplings[index] << std::endl;
            }
        }
        for(unsigned int i = this->num_fields; i < this->num_fields_and_couplings; ++i){
            fields_and_couplings[i] = 0.0;

        }
    
}


//Map pair-site frequencies indices
unsigned int PlmDCA::mapIndexPairSiteFreqs(const unsigned int i, const unsigned int j, 
    const unsigned int a, const unsigned int b)
{
    /*Maps pair site indices into contiguous memory layout

        Parameters
        ----------
            this    : PlmDCA instance
            i       : Trailing site in sequence
            j       : Following site in sequence
            a       : Residue at trailing site 
            b       : Residue at following site

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
            i       : Trailing index  when refering Jij(a,b) 
            j       : Following index when refering to Jij(a, b)
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
   
    unsigned int k = a + i * this->num_site_states; 
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
            sites

        a   : Residue at the trailing site 
        b   : Residue at the following site

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
        *instance   : PlmDCA class pointer
        x           :
        gradient    : 
        n           :
        step        :

    Returns
    -------
        void        : void

    */
   

    auto const& q = this->num_site_states;
    auto const& L = this->seqs_len;
    auto const& Nseq = this->num_seqs;
    auto const& lh = this->lambda_J;
    auto const& lJ = this->lambda_J;

    float fx = 0.0;
    for(unsigned int i = 0; i < this->num_fields_and_couplings; ++i) grad[i] =  0.0;


    // Compute gradients of the negative of pseudolikelihood from alignment data. 
    #pragma omp parallel for num_threads(this->num_threads)
    for(unsigned int i = 0; i < L; ++i){
        std::vector<float> prob_ni(q, 0.0);
        std::vector<float> fields_gradient(q, 0.0);
        std::vector<float> couplings_gradient(L * q * q, 0.0);
        float  fxi = 0.0;
        for(unsigned int n = 0; n < Nseq; ++n){
            // compute probability at site i for sequence n
            auto const& current_seq = this->seqs_int_form[n];
            for(unsigned int a = 0; a < q; ++a){
                auto const& indx_ia = this->mapIndexFields(i,a);
                prob_ni[a] += fields_and_couplings[indx_ia];
            }
            
            for(unsigned int j = 0; j < i; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    auto const& indx_jia = this->mapIndexCouplings(j, i, current_seq[j] - 1, a);
                    prob_ni[a] += fields_and_couplings[indx_jia];
                }
            }

            for(unsigned int j = i + 1; j < L; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    auto const&  indx_ija = this->mapIndexCouplings(i, j, a, current_seq[j] - 1);
                    prob_ni[a] += fields_and_couplings[indx_ija];
                }
            }
            auto max_exp = *std::max_element(prob_ni.begin(), prob_ni.end());
            for(unsigned int a = 0; a <  q; ++a) prob_ni[a] = std::exp(prob_ni[a] - max_exp);
            
            float zni = std::accumulate(prob_ni.begin(), prob_ni.end(), 0.0);
            zni = 1.0/zni;
            for(unsigned int a = 0; a < q; ++a) prob_ni[a] *= zni;
            
            //compute the gradient of the minus log likelihood with respect to fields and couplings
            auto current_seq_weight = this->seqs_weight[n];
            auto const& res_i = current_seq[i] - 1;

            fxi -= current_seq_weight * std::log(prob_ni[res_i]);

            fields_gradient[res_i] -= current_seq_weight;
            for(unsigned int a = 0; a < q; ++a) fields_gradient[a] += current_seq_weight * prob_ni[a];
            
            for(unsigned int j = 0; j < i; ++j){
                auto const& indx_ji = this->mapIndexCouplingsOneSite(j, current_seq[j] - 1, res_i);
                 
                couplings_gradient[indx_ji]  -= current_seq_weight;
            }
            for(unsigned int j = i + 1; j < L; ++j){
                auto const& indx_ij = this->mapIndexCouplingsOneSite(j, res_i,  current_seq[j] - 1);
                couplings_gradient[indx_ij] -= current_seq_weight;
            }

            for(unsigned int j = 0; j < i; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    auto const& indx_jia = this->mapIndexCouplingsOneSite(j, current_seq[j] -1, a);
                    couplings_gradient[indx_jia] += current_seq_weight * prob_ni[a];
                }
            }
            for(unsigned int j = i + 1; j < L; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    auto const& indx_ija = this->mapIndexCouplingsOneSite(j, a, current_seq[j] -1);
                    couplings_gradient[indx_ija] += current_seq_weight * prob_ni[a];
                }
            }   
        }// End of iteration through sequences
        
        #pragma omp critical
        {
            fx += fxi;
            for(unsigned int a = 0; a < q; ++a){
                auto const& indx_ia = this->mapIndexFields(i, a);
                grad[indx_ia] += fields_gradient[a];

            }

            for(unsigned int j = 0; j < i; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    for(unsigned int b = 0; b < q; ++b){
                        auto const& indx_jiab = this->mapIndexCouplings(j, i, a, b);
                        grad[indx_jiab] += couplings_gradient[this->mapIndexCouplingsOneSite(j, a, b)];
                    }
                }
            }
            for(unsigned int j = i + 1; j < L; ++j){
                for(unsigned int a = 0; a < q; ++a){
                    for(unsigned int b = 0; b < q; ++b){
                        auto const& indx_ijab = this->mapIndexCouplings(i, j, a, b);
                        grad[indx_ijab] += couplings_gradient[this->mapIndexCouplingsOneSite(j, a, b)];
                        
                    }
                }
            }
        }//end of critical block
    
    } // End of iteration through sites

    // Gradients of L2-norm regularization terms.
    for(unsigned int i = 0; i < L; ++i){
        for(unsigned int a = 0; a < q; ++a){
            auto const& indx_ia = this->mapIndexFields(i, a);
            auto const& hia = fields_and_couplings[indx_ia];
            grad[indx_ia] += 2.0 * lh * hia;
            fx += lh *  hia * hia;
        }
    }

    for(unsigned int i = 0; i < L - 1; ++i){
        for(unsigned int j = i + 1; j < L; ++j){
            for(unsigned int a = 0; a < q; ++a){
                for(unsigned int b = 0; b < q; ++b){
                    auto const& indx_ij_ab = this->mapIndexCouplings(i, j, a, b);
                    auto const& Jijab = fields_and_couplings[indx_ij_ab];
                    grad[indx_ij_ab] += 2.0 * lJ * Jijab;
                    fx += lJ *  Jijab *  Jijab;
                }
            }
        }
    }
    return fx;
 }


/*
void PlmDCA::printDCAScores(VectorXf const& h_and_J)
{   

    for(unsigned int i = 0; i < this->seqs_len - 1;  ++i){
        for(unsigned int j = i + 1; j < this->seqs_len; ++j){
            float scoreij = 0.0;
            for(unsigned int a = 0; a < this->num_site_states; ++a){
                for(unsigned int b = 0; b <  this->num_site_states; ++b){
                    if(a == this->num_site_states - 1 || b == this->num_site_states - 1)continue;
                    auto const& Jijab = h_and_J[this->mapIndexCouplings(i, j, a, b)];
                    scoreij += Jijab * Jijab;
                }
            }
            std::cout  << i + 1 << "\t" << j + 1 << "\t" << std::sqrt(scoreij) << std::endl; 
        }
    }
}
*/

//compute sequences weight
std::vector<float> PlmDCA::computeSeqsWeight()
{
    /*  Computes sequences weight.

        Parameters
        ----------
            this    : An instance of PlmDCA class.

        Returns
        -------
            m_seqs_weight   : The weigh of sequences computed with respect to 
                sequence identity obtained from this->seqid.

    */
    std::vector<float> m_seqs_weight(this->num_seqs);
    //std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    //std::cout << "#computing sequences weights" << std::endl;

    #if defined(_OPENMP)
        //initialize weights to zero for the parallel case since each sequences is going to be compared with itself.
        for(unsigned int i=0; i < this->num_seqs; ++i){
            m_seqs_weight[i] = 0.0; 
        }
        #pragma omp parallel for num_threads(this->num_threads)
        for(unsigned int i = 0; i < this->num_seqs; ++i){
            //seq_i = seqs_int_form[i];
            //printf("Thread: %d\n", omp_get_thread_num());
            float similarity_ij;
            unsigned int num_identical_residues;
            for(unsigned int j = 0; j < this->num_seqs; ++j){
                num_identical_residues = 0;
                for(unsigned int site =0; site < this->seqs_len; ++site){
                    if(this->seqs_int_form[i][site] == this->seqs_int_form[j][site]) num_identical_residues++;
                }
                similarity_ij = (float)num_identical_residues/(float)this->seqs_len;
                if (similarity_ij > this->seqid) m_seqs_weight[i] += 1.0;
            }
        }
    #else
        // Initialize weights to one for unique sequence pair comparisons.
        for(unsigned int i = 0; i < this->num_seqs; ++i){
            m_seqs_weight[i] = 1.0;
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
                    m_seqs_weight[i] += 1.0;
                    m_seqs_weight[j] += 1.0;
                }
            }
        }
    #endif
    
    //"Normalize" sequences weight 
    for(unsigned int i = 0; i < this->num_seqs; ++i) m_seqs_weight[i] = 1.0/m_seqs_weight[i];
    //for(unsigned int i = 0; i < this->num_seqs;++i) std::cout << i << " " << m_seqs_weight[i] << std::endl;
    //std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    return m_seqs_weight;
}



//print sequences for debugging 
void PlmDCA::printWeights()
{
    for(unsigned int i = 0; i < this->num_seqs; ++i){
        std::cout << "# " << i + 1 << " " << this->seqs_weight[i] << std::endl;
    }
    std::cout << "#Meff = " << std::accumulate(this->seqs_weight.begin(), this->seqs_weight.end(), 0.0) << std::endl;
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
    /*  Protein residues mapping
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W':19, 'Y':20,
        '-':21, '.':21, '~':21,
    */
        res_mapping['A'] = 1; res_mapping['C'] = 2; res_mapping['D'] = 3;
        res_mapping['E'] = 4; res_mapping['F'] = 5; res_mapping['G'] = 6;
        res_mapping['H'] = 7; res_mapping['I'] = 8; res_mapping['K'] = 9;
        res_mapping['L'] = 10; res_mapping['M'] = 11; res_mapping['N'] = 12;
        res_mapping['P'] = 13; res_mapping['Q'] = 14; res_mapping['R'] = 15;
        res_mapping['S'] = 16; res_mapping['T'] = 17; res_mapping['V'] = 18;
        res_mapping['W'] = 19; res_mapping['Y'] = 20; res_mapping['-'] = 21;
        res_mapping['.'] = 21; res_mapping['~'] = 21; res_mapping['B'] = 21;
        res_mapping['J'] = 21; res_mapping['O'] = 21; res_mapping['U'] = 21;
        res_mapping['X'] = 21; res_mapping['Z'] = 21;
    }else{
    /* RNA residues mapping
        'A':1, 'C':2, 'G':3, 'U':4, '-':5, '.':5, '~':5
    */
       res_mapping['A'] = 1; res_mapping['C'] = 2; res_mapping['G'] = 3;
       res_mapping['U'] = 4; res_mapping['-'] = 5; res_mapping['~'] = 5;
       res_mapping['.'] = 5; res_mapping['B'] = 5; res_mapping['D'] = 5;
       res_mapping['E'] = 5; res_mapping['F'] = 5; res_mapping['H'] = 5;
       res_mapping['I'] = 5; res_mapping['J'] = 5; res_mapping['K'] = 5;
       res_mapping['L'] = 5; res_mapping['M'] = 5; res_mapping['N'] = 5;
       res_mapping['O'] = 5; res_mapping['P'] = 5; res_mapping['Q'] = 5;
       res_mapping['R'] = 5; res_mapping['S'] = 5; res_mapping['T'] = 5;
       res_mapping['V'] = 5; res_mapping['W'] = 5; res_mapping['X'] = 5;
       res_mapping['Y'] = 5; res_mapping['Z'] = 5;
    }
    std::vector<std::vector<unsigned int>> seqs_int_form;
    std::ifstream msa_file_stream(this->msa_file);
    std::string current_line;
    std::string current_seq;
    int line_counter=0;
    int unique_seq_counter = 0;
    std::vector<unsigned int> current_seq_int;

    
    if(msa_file_stream.fail()){
        std::cerr << "Unable to open file " << this->msa_file << std::endl;
        exit(-1);
    }

    while(std::getline(msa_file_stream, current_line)){
        if(!current_line.empty() && current_line[0] != '>'){ 
            for(unsigned int i = 0; i < seqs_len; ++i){
                
                current_seq_int.emplace_back(res_mapping.at(toupper(current_line[i])));
            }
            //current_seq_int.shrink_to_fit();
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
    /*
    std::cout << "#" << __PRETTY_FUNCTION__ << std::endl;
    std::cout << "#Total number of sequences found: " << line_counter << std::endl;
    std::cout << "#Total number of unique sequences found: " << unique_seq_counter << std::endl;
    std::cout << "#" <<__PRETTY_FUNCTION__ << std::endl;
    */
    return seqs_int_form;
}
