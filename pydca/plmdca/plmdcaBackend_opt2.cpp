#include "plmdca_opt2.h"
#include "LBFGS.h"


using Eigen::VectorXf;
using namespace LBFGSpp;


extern "C" float* plmdcaBackend(unsigned short const biomolecule, unsigned short const num_site_states, 
    const char* msa_file, unsigned int const seqs_len, float const seqid, float const lambda_h, 
    float const lambda_J, unsigned int const max_iteration)
{  
    
    //std::cout << "#***********--PLMDCA BACKEND--**************" << std::endl;
    const int total_num_params = seqs_len * num_site_states + seqs_len * (seqs_len - 1) * num_site_states * num_site_states/2 ;
     
    static PlmDCA plmdca_inst(msa_file, biomolecule, seqs_len, num_site_states, seqid, lambda_h, lambda_J);

    class PlmDCALBFGSPP{
    private:
        int n;
    public: 
        PlmDCALBFGSPP(int m_n):n(m_n) {}
        float operator()(const VectorXf x, VectorXf& grad)
        {
            auto fx = plmdca_inst.gradient(x, grad);
            return fx;
        }

    };
    //Eigen::setNbThreads(1);
    
    const int N = total_num_params;
    // Set up parameters
    LBFGSParam<float> param;
    param.epsilon = 1e-6;
    param.max_iterations = max_iteration;
    // Create solver and function object
    LBFGSSolver<float> solver(param);
    PlmDCALBFGSPP fun(N);

    // Initial guess
    VectorXf x(N); // = VectorXf::Zero(n);
    plmdca_inst.initFieldsAndCouplings(x);
    // x will be overwritten to be the best point found
    float fx;
    int niter = solver.minimize(fun, x, fx);
    //plmdca_inst.printDCAScores(x);
    std::cerr << niter << " iterations" << std::endl;
    //std::cout << "x = \n" << x.transpose() << std::endl;
    std::cerr << " f(x) = " << fx << std::endl;
   
  float * h_and_J = nullptr;
    try{
        h_and_J = new float [total_num_params];
    }catch(std::bad_alloc &e){
        std::cerr << "Unable to allocate memory for fields and couplings backend\n";
        std::cerr << e.what() << "\n";
    }catch(std::exception &e){
        std::cerr << "Problem encountered in fields and couplings backend\n";
        std::cerr << e.what() << "\n";
    }
    for(int i = 0; i < total_num_params; ++i) h_and_J[i] = x[i];
    
    //std::cout << "#***********--PLMDCA BACKEND--**************" << std::endl;
    return h_and_J;
}


extern "C" void freeFieldsAndCouplings(void * h_and_J)
{  
    /*  Frees memory that has been used to store fields and couplings before 
        they are captured in the Python interface. 
        Parameter h_and_J must pass using ctypes.byref from Python.

    */
   float* h_and_J_casted = static_cast<float*>(h_and_J);  
    if(h_and_J_casted !=nullptr){
        delete [] h_and_J_casted;
        h_and_J_casted = nullptr;
        //std::cerr << __PRETTY_FUNCTION__ << "\n";
        //std::cerr<< "DELETED fields_and_couplings from the C++ Backend\n";
    }
}
    