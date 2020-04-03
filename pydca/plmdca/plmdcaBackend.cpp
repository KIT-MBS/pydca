#include "include/plmdca.h"

/*Implements the pseudolikelihood maximization direct couplings analysis 
for protein and RNA sequences.

Authors: Mehari B. Zerihun, Fabrizio Pucci

*/


class ObjectiveFunction{
    /*Objective Function for lbfgs input. 
        
    Attributes
    ----------
        m_x     : A dynamic array containing fields and couplings.
        plmdca_inst : PlmDCA 
        m_verbose : bool 
        m_max_iterations : unsigned int 
    */
   
    protected:
        float* m_x;
        bool m_verbose;
        unsigned int m_max_iterations;
        PlmDCA plmdca_inst;
        
    public:
        ObjectiveFunction(unsigned short const biomolecule, unsigned short const num_site_states, 
            const char* msa_file, unsigned int const seqs_len, float const seqid, float const lambda_h, 
            float const lambda_J, unsigned int const max_iterations, const unsigned int num_threads, bool verbose):
            m_x(NULL), 
            m_verbose(verbose), 
            m_max_iterations(max_iterations), 
            plmdca_inst(msa_file, biomolecule, seqs_len, num_site_states, seqid, lambda_h, lambda_J, num_threads)
        {
            // ObjectiveFunction constructor body
        }        


        float* getFieldsAndCouplings() 
        {
            return this->m_x; 
        }


        int run(int N)
        {
            /*Performs plmDCA computation using LBFGS optimization.

            Parameters
            ----------
                N       : Total number of fields and couplings. 

            Returns
            -------
                ret     : Exit status of LBFGS optimization.
            */
        
            float fx;
            this->m_x = lbfgs_malloc(N);

            if (this->m_x == NULL) {
                printf("ERROR: Failed to allocate a memory block for variables.\n");
                return 1;
            }
            //initialize parameters
            lbfgs_parameter_t param;
            lbfgs_parameter_init(&param);
            param.epsilon = 1E-3;
            param.max_iterations = this->m_max_iterations;
            param.max_linesearch = 5;
            param.ftol = 1E-4;
            //param.wolfe = 0.2;
            param.m = 5 ;

            this->plmdca_inst.initFieldsAndCouplings(m_x);
            //Start the L-BFGS optimization; this will invoke the callback functions
            //evaluate() and progress() when necessary.
                        
            int ret = lbfgs(N, m_x, &fx, _evaluate, _progress, this, &param);
            // return status value ret == -1001 corresponds with convergence for a given precision
            /* Report the result. */
            if(this->m_verbose){
                if (ret==-1001){
                    fprintf(stderr, "L-BFGS optimization completed\n");
                }else{
                    fprintf(stderr, "L-BFGS optimization terminated with status code = %d\n", ret);
                    fprintf(stderr, "fx = %f\n", fx);
                }
            }
        
            return ret;
        }


    protected:
        static float _evaluate( void* instance, const float*x, float* g, const int n, const float step)
        {
            /*Computes the gradient of the regularized negative pseudolikelihood function for 
            protein/RNA alignments. 

            Parameters
            ----------
                instance    : An instance of ObjectiveFunction class. 
                x           : Array of fields and couplings.
                g           : Array of gradient of the negative log pseudolikelihood
                    of the conditional probablity for protein/RNA alignments.
                n           : Number of fields and couplings?
                step        : The step size for gradient decent.

            Returns
            --------
                fx          : Value of plmDCA objective function

            */

            return reinterpret_cast<ObjectiveFunction*>(instance)->evaluate(x, g, n, step);
        }


        float evaluate(const float*x, float* g, const int n, const float step)
        {
            float fx;
            fx = this->plmdca_inst.gradient(x, g);
            return fx;
        }


        static int _progress(void* instance, const float* x, const float* g, const float fx, 
            const float xnorm, const float gnorm, const float step, int n, int k, int ls)
        {
            return reinterpret_cast<ObjectiveFunction*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
        }


        int progress(const float* x, const float* g, const float fx, const float xnorm, const float gnorm,
            const float step, int n, int k, int ls)
        {
            if(this->m_verbose){
                fprintf(stderr, "Iteration %d:\n", k);
                fprintf(stderr, "fx = %f, xnorm = %f, gnorm = %f, step = %f\n", fx, xnorm, gnorm, step);
                fprintf(stderr, "\n");
            }
            return 0;
        }
};
    


extern "C" float* plmdcaBackend(unsigned short const biomolecule, 
    unsigned short const num_site_states, 
    const char* msa_file, unsigned int const seqs_len, 
    float const seqid, float const lambda_h, 
    float const lambda_J, unsigned int const max_iteration, 
    const unsigned int num_threads, bool verbose )
{  
    /*Interface for the Python implementation of plmDCA. 

    Parameters
    ----------
        biomolecule     : Type of biomolecule (protein or RNA).
        num_site_states : Number of states/residues plus gap for MSA data.
        msa_file        : Path to the FASTA formatted MSA file.
        seqs_len        : The length of sequences in MSA data.
        seqid           : Sequence identity threshold.
        lambda_h        : Regularization parameter for fields.
        lambda_J        : Regularization parameter for couplings.
        max_iteration   : Maximum number of gradient decent iterations.
        num_threads     : Number of threads for PlmDCA (when OpenMP is supported).
        verbose         : Print logging message on the terminal.

    Returns
    -------
        h_and_J        : Fields and couplings array. This data is fetched into the
            Python interface. 

    */
   #if defined(_OPENMP)
    // can use multiple threads
    #else 
        if(num_threads > 1){
            std::cerr << "Cannot set multiple threads when OpenMP is not supported\n";
            throw std::runtime_error("Invalid number of threads");
        }
    #endif

    const int total_num_params = seqs_len * num_site_states + seqs_len * (seqs_len - 1) * num_site_states * num_site_states/2 ; 

    // Start computation 
    ObjectiveFunction objfun_inst(biomolecule, num_site_states, msa_file, seqs_len, seqid, 
        lambda_h, lambda_J, max_iteration, num_threads, verbose
    );

    //const int N = total_num_params;

    objfun_inst.run(total_num_params);
    auto h_and_J = objfun_inst.getFieldsAndCouplings();

    return h_and_J;
}


extern "C" void freeFieldsAndCouplings(void * h_and_J)
{  
    /*  Frees memory that has been used to store fields and couplings before 
        they are captured in the Python interface.

        Parameters
        ----------
            h_and_J : Pointer to the fields and couplings vector 
        
        Returns
        -------
            void    : No return value

    */
   float* h_and_J_casted = static_cast<float*>(h_and_J);  
    if(h_and_J_casted !=nullptr){
        delete [] h_and_J_casted;
        h_and_J_casted = nullptr;
    }
}
    