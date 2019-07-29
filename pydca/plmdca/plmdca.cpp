#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <locale>
#include <omp.h>

/*	Implementation of pseudo-likelihood maximization direct coupling analysis (plmDCA) 
 *	
 *	Author Emanuel Karl Peter 
 *
 */

int main(int argc, char **argv) {

    
    int n,i,j,k,l;
    
    if (argc < 2) {
        std::cerr << " Wrong format: " << argv[0] << " [infile] " << std::endl;
        return -1;
    }

    std::ifstream input(argv[1]);
    if (!input.good()) {
        std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
        return -1;
    }
    std::string line, id, DNA_sequence;
    std::string String;
    std::vector <std::string > Sequence_vector;
    char  inputstring[100];
    double tolerance=1E-6;
    int  **msa_vector,length;
    double ***msa_vec_d,***omega_d;
    double **single_frequency_count;
    double ***double_frequency_count;
  //  double double_frequency_count[150][150][6][6]; 
    double *weight;
    double seq_id=0.8;
    double q_parameter;
    double pseudocount = 0.5;
    double *G_of_r,**h,****J;
    double **logPot;
    int    n2;
    int    maxcontact = 1000;
    q_parameter = 5.0;
    double norm;
    int number_of_min_steps = 200;
    int num_iterations = 2;
    double **d_G_dh,****d_G_dJ;
    int    switchflag=1;
    double stepsize_H = 1E+1;
    double stepsize_J = 1E+1;
    double threshold_B=0.95;
      
    int q_int = 5;
    
    n = 0;
    
     if(argc < 7 && argc > 2) {
        
        std::cout << "Program num_min_steps num_iter tol maxcontact [rna or protein]" << "\n";

        return 0;
        
    };   
    
#if defined(_OPENMP)    
#pragma omp parallel 
#endif       
    
    
#pragma omp critical    
    if(argc > 6) {
    
    strcpy(inputstring,argv[2]);
    sscanf(inputstring,"%d",&number_of_min_steps);
    strcpy(inputstring,argv[3]);
    sscanf(inputstring,"%d",&num_iterations);
    strcpy(inputstring,argv[4]);
    sscanf(inputstring,"%lg",&tolerance);  
    strcpy(inputstring,argv[5]);
    sscanf(inputstring,"%d",&maxcontact);   
    strcpy(inputstring,argv[6]);    
    if(strcmp(inputstring,"rna") == 0) switchflag = 1;
    if(strcmp(inputstring,"protein") == 0) switchflag = 2;
    
   };

   if(switchflag == 2) {
       q_int = 21;
       q_parameter = 21.0;  
   };   
    
    // Parse MSA, argv[1]
    
      while (std::getline(input, line)) {

            id = line.substr(1);
            DNA_sequence.clear();
            DNA_sequence += line;     
            
         if(line[0] != '>') { 
            if(!id.empty())
                std::cout << id << " : " << DNA_sequence << std::endl;
            if(!id.empty()) {
                std::cout << DNA_sequence.length() << " seq. length " << std::endl;
                Sequence_vector.push_back(DNA_sequence);
                length = DNA_sequence.length();
                n++;
            };
         }; 
      };
      
    //  return 0;
      
    // Initialize MSA Vector
    
      std::cout << n << " Number of lines " << std::endl;
      
      int number_of_sequences = n-1;
      
      n = 0;
      while(n <= number_of_sequences) {
        
          std::cout << Sequence_vector[n] << std::endl;
          n++;
          
      };
      
      std::cout << length << " length  " << number_of_sequences << " number of sequences "<< std::endl;
      
      msa_vector = (int **) malloc(sizeof(int)*(number_of_sequences*2));
      
      for(n=0;n<=number_of_sequences+1;n++) {
        
          msa_vector[n] = (int *) malloc(sizeof(int)*(length+1));
          
      };
      
      msa_vec_d = (double ***) malloc(sizeof(double)*(number_of_sequences*2));
      
      for(n=0;n<=number_of_sequences+1;n++) {
        
          msa_vec_d[n] = (double **) malloc(sizeof(double)*(length+1));
          
          for(n2=0;n2<=length;n2++) {
            
              msa_vec_d[n][n2] = (double *) malloc(sizeof(double)*(q_int+1));
              
          };
      };   
      
      omega_d = (double ***) malloc(sizeof(double)*(number_of_sequences*2));
      
      for(n=0;n<=number_of_sequences+1;n++) {
        
          omega_d[n] = (double **) malloc(sizeof(double)*(length+1));
          
          for(n2=0;n2<=length;n2++) {
            
              omega_d[n][n2] = (double *) malloc(sizeof(double)*(q_int+1));
              
          };
      };        
      
      for(n=0;n<=number_of_sequences*2-1;n++) {
          for(i=0;i<String.length();i++) {
            
              msa_vector[n][i] = 0;
              
              for(n2=0;n2<=q_int;n2++) {
                
                  msa_vec_d[n][i][n2] = 0.0;
                  
              };
          };
      };
      
     // Set nucleotide indices for each entry in each sequence in a matrix Nxi
      
      n = 0;
      
      while(n <= number_of_sequences) {
          
          String.clear();
          String = Sequence_vector[n];
          
          for(i=0;i<String.length();i++) {
              
         //     std::cout << String.length() << " i " << i << " n " << n << " str.len " << std::endl;
              
            if(switchflag == 1) {  
              
              if(String[i] == 'A') {
                  msa_vector[n][i] = 1;
                  msa_vec_d[n][i][1] = 1.0;
              };
              if(String[i] == 'G') {
                  msa_vector[n][i] = 2;
                  msa_vec_d[n][i][2] = 1.0;                  
              };
              if(String[i] == 'C') {
                  msa_vector[n][i] = 3;
                  msa_vec_d[n][i][3] = 1.0;                  
              };
              if(String[i] == 'U') {
                  msa_vector[n][i] = 4;
                  msa_vec_d[n][i][4] = 1.0;                  
              };
              if(String[i] == '-') {
                  msa_vector[n][i] = 5; 
                  msa_vec_d[n][i][5] = 1.0;                  
              };
         //     std::cout << msa_vector[n+1][i+1];
            };
            
             if(switchflag == 2) {  
              
              if(String[i] == 'A') {
                  
                  msa_vector[n][i] = 1;
                  msa_vec_d[n][i][1] = 1.0;
              };
              if(String[i] == 'C') {
                  msa_vector[n][i] = 2;
                  msa_vec_d[n][i][2] = 1.0;
              }
              if(String[i] == 'D') {
                  msa_vector[n][i] = 3;
                  msa_vec_d[n][i][3] = 1.0;
              }    
              if(String[i] == 'E') {
                  msa_vector[n][i] = 4;
                  msa_vec_d[n][i][4] = 1.0;
              }    
              if(String[i] == 'F') {
                  msa_vector[n][i] = 5; 
                  msa_vec_d[n][i][5] = 1.0;
              }    
              if(String[i] == 'G') {
                  msa_vector[n][i] = 6;
                  msa_vec_d[n][i][6] = 1.0;
              }   
              if(String[i] == 'H') {
                  msa_vector[n][i] = 7;
                  msa_vec_d[n][i][7] = 1.0;
              }    
              if(String[i] == 'I') {
                  msa_vector[n][i] = 8;
                  msa_vec_d[n][i][8] = 1.0;
              }   
              if(String[i] == 'K') {
                  msa_vector[n][i] = 9;
                  msa_vec_d[n][i][9] = 1.0;
              }    
              if(String[i] == 'L') {
                  msa_vector[n][i] = 10; 
                  msa_vec_d[n][i][10] = 1.0;
              }    
              if(String[i] == 'M') {
                  msa_vector[n][i] = 11;
                  msa_vec_d[n][i][11] = 1.0;
              }    
              if(String[i] == 'N') {
                  msa_vector[n][i] = 12;
                  msa_vec_d[n][i][12] = 1.0;
              }    
              if(String[i] == 'P') {
                  msa_vector[n][i] = 13;
                  msa_vec_d[n][i][13] = 1.0;
              }    
              if(String[i] == 'Q') {
                  msa_vector[n][i] = 14;
                  msa_vec_d[n][i][14] = 1.0;
              }   
              if(String[i] == 'R') {
                  msa_vector[n][i] = 15; 
                  msa_vec_d[n][i][15] = 1.0;                  
              }    
              if(String[i] == 'S') {
                  msa_vector[n][i] = 16;
                  msa_vec_d[n][i][16] = 1.0;                  
              }    
              if(String[i] == 'T') {
                  msa_vector[n][i] = 17;
                  msa_vec_d[n][i][17] = 1.0;
              }    
              if(String[i] == 'V') {
                  msa_vector[n][i] = 18;
                  msa_vec_d[n][i][18] = 1.0;
              }    
              if(String[i] == 'W') {
                  msa_vector[n][i] = 19;
                  msa_vec_d[n][i][19] = 1.0;
              }    
              if(String[i] == 'Y') {
                  msa_vector[n][i] = 20; 
                  msa_vec_d[n][i][20] = 1.0;
              }    
              if(String[i] == '-') {
                  msa_vector[n][i] = 21; 
                  msa_vec_d[n][i][21] = 1.0;
               }   
         //     std::cout << msa_vector[n+1][i+1];
         
              for(k=1;k<=q_int;k++) {
                 
                  omega_d[n][i][k] = 0.0;
                  
               };
         
            };           
              
          };
          
        //  std::cout << "\n";
          
          n++;
      };
      
      weight = (double *) malloc(sizeof(double)*number_of_sequences+1);
      
      // Calculate sequence weights 
      
      double id2;
      
      std::cout << " Calculate sequence weights over " << number_of_sequences << " sequences. Threshold :: "<< seq_id << std::endl;
      
      for(i=0;i<=number_of_sequences;i++) {
        
          weight[i] = 0.0;
          
      };
      
      for(i=0;i<=number_of_sequences;i++) {
          
          weight[i] += 1.0;
          
          for(j=0;j<=number_of_sequences;j++) {
             
            if(i != j) {  
                
              id2 = 0.0;  
              
              for(k=0;k<=length;k++) {
                  
               if(msa_vector[i][k] == msa_vector[j][k]) id2 = id2 + 1.0;       
                  
               };
           };
          
          id2 = id2/((double)length);
          
       //   std::cout << " id : " << id2 << " index " << i << std::endl; 
          weight[i] += 1.0/id2;
          weight[j] += 1.0/id2;
          
          if(id2 >= (1.0 - seq_id)) {
            
              weight[i] += 1.0;
              weight[j] += 1.0;
              
          };
        };
      };
      
    //  Sequence_vector.clear();
     
      double B_eff,mean,max = 0.0,min=1000.0;
      
      B_eff = 0.0;
      mean  = 0.0;
      
      for(i=0;i<=number_of_sequences;i++) {
          
          weight[i] = 1.0/weight[i];
          B_eff += weight[i];
                  
        //  std::cout << weight[i] << " " << mean << std::endl;
          
      };
      
      double lambda_J,lambda_h;
      double scaled_lambda_h,scaled_lambda_J;
      
    if(B_eff>500.0) {
        
        lambda_J=0.01;
        
    } else {
        
        lambda_J=0.1-(0.1-0.01)*B_eff/500.0;
        
    };
    
    lambda_h = 0.01;
    
    lambda_h=lambda_J;
    scaled_lambda_h=lambda_h*B_eff;   
    scaled_lambda_J=lambda_J*B_eff/2.0;      

      std::cout << " Generation of single and double site couplings h_i, J_ij " << std::endl;
      
      h      = (double **) malloc( sizeof(double) * (length+1));
      
      for(k=0;k<=length;k++) {
        
          h[k] = (double *) malloc( sizeof(double) * (q_int+1));
          
      };
      
      d_G_dh      = (double **) malloc( sizeof(double) * (length+1));
      
      for(k=0;k<=length;k++) {
        
          d_G_dh[k] = (double *) malloc( sizeof(double) * (q_int+1));
          
      };      
      
      J = (double ****) malloc( sizeof(double) * (length+1));
          
      for(k=0;k<=length;k++) {
              
         J[k] = (double ***) malloc( sizeof(double) * (length+1));
              
         for(n=0;n<=length;n++) {
          
             J[k][n] = (double **) malloc( sizeof(double) * (q_int+1));
             
          for(n2=0;n2<=q_int;n2++) {
            
              J[k][n][n2] = (double *) malloc( sizeof(double) * (q_int + 1)); 
              
          };
        };
      };  
      
    d_G_dJ = (double ****) malloc( sizeof(double) * (length+1));
          
      for(k=0;k<=length;k++) {
              
         d_G_dJ[k] = (double ***) malloc( sizeof(double) * (length+1));
              
         for(n=0;n<=length;n++) {
          
             d_G_dJ[k][n] = (double **) malloc( sizeof(double) * (q_int+1));
             
          for(n2=0;n2<=q_int;n2++) {
            
              d_G_dJ[k][n][n2] = (double *) malloc( sizeof(double) * (q_int + 1)); 
              
          };
        };
      };       
      
   std::cout << " here " << std::endl;
         
      for(k=0;k<=length;k++) {
          
        for(n=1;n<=q_int;n++) {  
          
             h[k][n] = 0.0;
        
        for(l=0;l<=length;l++) {              
          
            for(n2=1;n2<=q_int;n2++) {
            
              J[k][l][n][n2] = 0.0;
        
            };
           };
         };
        };
     
  int    r;
  double **probability;
     
  probability = (double **) malloc(sizeof(double)* (length + 1));
      
  for(i=0;i<=length;i++) {
     
       probability[i] = (double *) malloc(sizeof(double)* (q_int+1));
       
   };
   
   double *Z,*G;
   
   Z = (double *) malloc(sizeof(double)* (length + 1)); 
   
   G = (double *) malloc(sizeof(double)* (length + 1));
   
   double sum_hJ,sum_GK,dG=0.0,dG2=0.0;
   int nn;   
   
   std::cout << " here " << std::endl;
   
 for(r=1;r<=number_of_min_steps;r++) {
     
  dG2 = 0.0;   
     
 for(nn = 1;nn <= num_iterations; nn++) {   
     
  for(k=0;k<=length;k++) {
      
    Z[k] = 0.0;  
                
    for(n=1;n<=q_int;n++) {   
          
        probability[k][n] = h[k][n];
        
        sum_hJ            = h[k][n];
         
        for(l=0;l<=length;l++) {                
                
           if(k < l) { 
             
             for(n2=1;n2<=q_int;n2++) {  
               
               probability[k][n] += J[k][l][n][n2];
               sum_hJ            += J[k][l][n][n2];
               
             };
           };
          };
      
        probability[k][n] = exp(probability[k][n]);
        Z[k]             += exp(sum_hJ);
          
      };
     };

   sum_GK = 0.0;  
     
   for(k=0;k<=length;k++) {
       
       G[k] = log(Z[k]);
       
       for(n=1;n<=q_int;n++) {  
           
          probability[k][n] /= Z[k];
          
          G[k] -= h[k][n];
          
        for(l=0;l<=length;l++) {                
                
           if(k < l) { 
             
             for(n2=1;n2<=q_int;n2++) { 
               
                 G[k] -= J[k][l][n][n2];
                 
             };
           };
          };
       //   std::cout << probability[k][n] << " prob " << "\n";          
          
      };
      
      sum_GK += G[k];
      
    };
   
   dG = 0.0;

   
for(i=0;i<=number_of_sequences;i++) {      
  //    if(weight[i] >= mean) {      

#if defined(_OPENMP)     
#pragma omp parallel for
#endif     
    
      for(k=0;k<=length;k++) {
          
         // lambda_h = 1.0/Z[k];
           
        for(n=1;n<=q_int;n++) {  
          
            d_G_dh[k][n] = -1.0/B_eff*weight[i]*((omega_d[i][k][n] + msa_vec_d[i][k][n]) - probability[k][n] + 2.0*(lambda_h)*h[k][n]);
            
          //  dG += d_G_dh;
          
            h[k][n] += d_G_dh[k][n]*stepsize_H;
            omega_d[i][k][n] -= d_G_dh[k][n]*probability[k][n]*1E-4;
            
         };  
   //  };  
   };
   };
   
   for(i=0;i<=number_of_sequences;i++) {
       
  //    if(weight[i] >= mean) {      
       
      for(k=0;k<=length;k++) {
          
         // lambda_h = 1.0/Z[k];
           
        for(n=1;n<=q_int;n++) { 
            
            dG += d_G_dh[k][n];
        }
      }
   }
   
   dG2 = dG/((double)length);
   
   std::cout << "Step h : " << nn << " DG ::" << dG << " GK :: " << sum_GK << "\n";
   
 };
        
 for(nn = 1;nn <= num_iterations; nn++) {   
     
  for(k=0;k<=length;k++) {
      
    Z[k] = 0.0;  
                
    for(n=1;n<=q_int;n++) {   
          
        probability[k][n] = h[k][n];
        
        sum_hJ            = h[k][n];
         
        for(l=0;l<=length;l++) {                
                
           if(k < l) { 
             
             for(n2=1;n2<=q_int;n2++) {  
               
               probability[k][n] += J[k][l][n][n2];
               sum_hJ            += J[k][l][n][n2];
               
             };
           };
          };
      
        probability[k][n] = exp(probability[k][n]);
        Z[k]             += exp(sum_hJ);
          
      };
     };

   sum_GK = 0.0;  
     
   for(k=0;k<=length;k++) {
       
       G[k] = log(Z[k]);
       
       for(n=1;n<=q_int;n++) {  
           
          probability[k][n] /= Z[k];
          
          G[k] -= h[k][n];
          
        for(l=0;l<=length;l++) {                
                
           if(k < l) { 
             
             for(n2=1;n2<=q_int;n2++) { 
               
                 G[k] -= J[k][l][n][n2];
                 
             };
           };
          };
       //   std::cout << probability[k][n] << " prob " << "\n";          
          
      };
      
      sum_GK += G[k];
      
    };
    
    dG = 0.0;

    
  for(i=0;i<=number_of_sequences;i++) {       
  //   if(weight[i] >= mean) {    
       
#if defined(_OPENMP)     
#pragma omp parallel for
#endif       
      
      for(k=0;k<=length;k++) {
          
        //  lambda_J = 1.0/Z[k];           
        for(n=1;n<=q_int;n++) {    
            
          if(msa_vec_d[i][k][n] == 1.0) { 
        
            for(l=0;l<=length;l++) {  
              
               for(n2=1;n2<=q_int;n2++) {
                
                d_G_dJ[k][l][n][n2] = -1.0/B_eff*weight[i]*(msa_vec_d[i][k][n])*(omega_d[i][l][n2] + msa_vec_d[i][l][n2] - probability[l][n2] + 2.0*(lambda_J)*J[k][l][n][n2]);
            
                J[k][l][n][n2] += d_G_dJ[k][l][n][n2]*stepsize_J;          
                
                omega_d[i][l][n2] -= d_G_dJ[k][l][n][n2]*probability[l][n2]*1E-4;
                
               };
              };
            };
           }; 
        };
   //    };
     };
     
      for(k=0;k<=length;k++) {
          
        //  lambda_J = 1.0/Z[k];           
        for(n=1;n<=q_int;n++) {    
        
            for(l=0;l<=length;l++) {  
              
               for(n2=1;n2<=q_int;n2++) {
                   
                   dG += d_G_dJ[k][l][n][n2];
                   
               };
              };
             };
            };
       
       if(nn == num_iterations) dG2 += dG/((double)(length*length));
       
       std::cout << "Step J : " << nn << " DG ::" << dG << " GK :: " << sum_GK << "\n";
       
      };
      
  dG2 /= (2.0*(double)number_of_sequences);    
   
  std::cout << "Step : " << r << " DG ::" << dG2 << " GK :: " << sum_GK << " tolerance " << tolerance << " number_of_min_steps " << number_of_min_steps << "\n";
       
  if(sqrt(pow(dG2,2)) <= tolerance) break;
  
  
  };
  
  std::cout << "Average over sequences " << std::endl;
    
     double **J_coupling;
     
     J_coupling = (double **) malloc(sizeof(double)*(length+1));
     
     for(k=1;k<=length;k++) {
        
         J_coupling[k] = (double *) malloc(sizeof(double)*(length+1));
         
      };
     
    double ***av_k,***av_l;
    double **av_kl;
    
    av_k = (double ***) malloc(sizeof(double)*(length+1));
    
     for(k=1;k<=length;k++) {
        
         av_k[k] = (double **) malloc(sizeof(double)*(length+1));
         
         for(l=1;l<=length;l++) {
           
             av_k[k][l] = (double *) malloc(sizeof(double)*(q_int+1));
             
         };
      };     
    
    av_l = (double ***) malloc(sizeof(double)*(length+1));
    
     for(k=1;k<=length;k++) {
        
         av_l[k] = (double **) malloc(sizeof(double)*(length+1));
         
         for(l=1;l<=length;l++) {
           
             av_l[k][l] = (double *) malloc(sizeof(double)*(q_int+1));
             
         };
      };    
    
    av_kl = (double **) malloc(sizeof(double)*(length+1));
     
     for(k=1;k<=length;k++) {
        
         av_kl[k] = (double *) malloc(sizeof(double)*(length+1));
         
      };    
    
    int n1 = 0;
    n2 = 0;
    
   for(k=1;k<=length;k++) {
        
     for(n=1;n<=q_int;n++) {      
        
      for(l=1;l<=length;l++) {
   
        for(n2=1;n2<=q_int;n2++) {    
            
          av_k[k][l][n2]  += J[k][l][n][n2]/((double)q_int);
          av_l[k][l][n]   += J[k][l][n][n2]/((double)q_int);
          
          av_kl[k][l]    += J[k][l][n][n2]/((double)(q_int*q_int));
          
        };
       };
     };
   };
    
    for(k=1;k<=length;k++) {
        
      for(n=1;n<=q_int;n++) {       
        
        for(l=1;l<=length;l++) {
   
          for(n2=1;n2<=q_int;n2++) {              
            
            J[k][l][n][n2] = J[k][l][n][n2] - av_l[k][l][n] - av_k[k][l][n2] + av_kl[k][l];
          
            J_coupling[k][l] += sqrt(pow(J[k][l][n][n2],2));
            
            };
          };
        };
    };
    
    double *av_k2,*av_l2;
    double av_kl2;   
    
    av_k2 = (double *) malloc(sizeof(double)*(length+1));
    av_l2 = (double *) malloc(sizeof(double)*(length+1));    
    
    for(k=1;k<=length;k++) {
            
      for(l=1;l<=length;l++) {
            
          av_k2[k]  += J_coupling[k][l]/((double)length);
          av_l2[l]  += J_coupling[k][l]/((double)length);
          
          av_kl2    += J_coupling[k][l]/((double)(length*length));
          
        };
       };   
    
    for(k=1;k<=length;k++) {
            
      for(l=1;l<=length;l++) {
      
          J_coupling[k][l] = J_coupling[k][l] - (av_k2[k]*av_l2[l])/av_kl2;
          
      };
    };
   
    for(k=1;k<=length;k++) {
            
      for(l=1;l<=length;l++) {
      
          J_coupling[k][l] = 0.5*(J_coupling[k][l] + J_coupling[l][k]);
          
      };
    };    
    
    int *k_index;
    int *l_index;
    double *J_index;
    
    r = 0;
    
    k_index = (int *) malloc(sizeof(int)*((length+1)*(length+1)));
    l_index = (int *) malloc(sizeof(int)*((length+1)*(length+1)));    
    J_index = (double *) malloc(sizeof(double)*((length+1)*(length+1)));     
    
    double buff_J1,buff_J2;
    int    buff_k1,buff_k2;
    int    buff_l1,buff_l2;
    
    for(k=1;k<=length;k++) {
        
        for(l=k;l<=length;l++) {
          
            k_index[r] = k;
            l_index[r] = l;
            J_index[r] = J_coupling[k][l];
         
         //   std::cout << length << " " << k << " " << l << " " << r << " r " << std::endl;
            
            r++;
            
        };
    };
    
    for(k=0;k<=r;k++) {
        for(l=0;l<=r;l++) 
        {
            if(J_index[l] < J_index[k]) {
               
                buff_J1 = J_index[k];
                buff_J2 = J_index[l];
                J_index[l] = buff_J1;
                J_index[k] = buff_J2;
                
                buff_k1  = k_index[k];
                buff_k2  = k_index[l];
                k_index[l] = buff_k1;
                k_index[k] = buff_k2;
                
                buff_l1  = l_index[k];
                buff_l2  = l_index[l];
                l_index[l] = buff_l1;
                l_index[k] = buff_l2;                
            }; 
        };
    };
        
    std::ofstream outfile,contact_file;
    
    outfile.open("DI_ij.dat");
    contact_file.open("contact_ij.dat");
  
    n = 1;
  
    FILE *fp = fopen("DI_ij.dat","w");    
    
     for(k=0;k<=r;k++) {
     
       if(k_index[k] != l_index[k]) {  
         
        fprintf(fp,"%d\t%d\t%15.25f\n",k_index[k],l_index[k],J_index[k]);       
         
        if(n <= maxcontact) {
        
            contact_file << std::setw(5) << k_index[k]+1 << " " << std::setw(5) << l_index[k]+1 << "\n";
        //    contact_file << std::setw(5) << l_index[k]+1 << " " << std::setw(5) << k_index[k]+1 << "\n";           
            
            n++;
        };
       }; 
     };
    
    fclose(fp); 
    contact_file.close();
    
    std::cout << " finalized " << std::endl;
      
}
