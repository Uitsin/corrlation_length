// From DL0115aseptest; the good1
// From From0120ModifyF
// From 0130Model1
// From Mf_noTDIP.cpp.cpp
// From newKexpr0222

// Test: No initial stretch. Will it break?  [Off]
// Note: Pure Neo-Hookean now
// X=100, Lf=7, ER=4. It's for know if the Lf should be propotional to a.
// [back to just ux] ux -> 4ux, ER and EL's conditions become 4*a
// //\\ Change energy functional into reciprocal gamma. (1113)
// //\\Modify the function "FuComponent" to let zero force is at uij = a. (1109) [not open at all]
// Xdim=Ydim=0.01
//// Add Lf0 and g for changing Lf.
// damping changed and time increments changed
// uijnorm output changed
// // Add I4's terms.
// K value of a complete hexagonal unit is 29.0725
// change the parameters K0 and g for Lf 
// Tried to make the incomplete hexagons subject to the same criterion.
// Modify M function with adding a parameter "onset" for the tubes to play a role.
// Still has the separate settings about GoUValue.
// From recovery0107
// T<=
// Delete the T=0's dissipation. 
// C3p4_o4p4



#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

using namespace std;

const double density = 930;
const double thick = 1.5e-4;
const double Xdim = 0.051;   // just arbitrary
const double Ydim = 0.051;

const int X = 172;   // size of (0.5, 0.866) direction
const int Y = 172;  // size of (0.5, -0.866) direction
const int NP = static_cast<int>( X%2==0 ? ((X+1)*(2*floor(0.5*X/1.7320508)+1)+2*X*(floor(0.5*X/1.7320508)+1)) : (X*(2*floor(0.5*X/1.7320508)+1)+2*(X+1)*floor(0.5*X/1.7320508)) );                           // Number of particles
const int T = 500000;
const double dt = 4e-7;
const double a = Xdim/X;                // original lattice constant. 
const double m = density*thick*Xdim*Ydim/NP;     
double Lf0 = 9;
double g = 600;
double K0 = 29.0725;
double A = 501;
double B = 0.053;
double beta = 9e-6;            // Lf is the lambda_f (>1)
const int ER = 24;      // Right Edge
const int EL = -24;     // Left Edge
const int SLL = static_cast<int>(ER - EL + 1);   // Slit Length
double Quotient = (2*floor(0.5*X/1.7320508)+1);
double stchX = 1;
double stchY = 5.45;
const double Axi = 0;                            // initial Ax and Ay
const double Ayi = 1;
double MAG = stchX * stchY;
const double CRIT = 3.8;    // Critical stretch
const double onset = 5.45;
const int tork = 4282;
int shu;   // last number of neighbor.
int time_no_snap = 100; 



double ENSQ (double x1, double y1, double x2, double y2) {     // Euclidean Norm Square.  Note: (x, y) (x, y)
    
    return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);

}
/*
double TDIP (double a1, double a2, double b1, double b2) {               // Arg: ( (a1, a2), (b1, b2) )
    
    return  a1*b1 + a2*b2;
    
}
*/
string int_to_string(int i) {       // Input an integer number and convert it into a sting.
  ostringstream sout;
  sout << i;
  return sout.str();
}

void UpdateList (int** nnl, double** uijx, double** uijy, double** uijNSQ, double** Ouijx, double** Ouijy, int* fn, int NP, int k, int j) {     // (nnl, uijx, uijy, uijNSQ, aijASQ, &fn, NP, k, j)
    
    int bk = *( *( nnl + k ) + j );
    int bj = -1;
              
    for (int s = j; s < (fn[k]-1); s++) {
             
         *( *(nnl + k) + s ) = *( *(nnl + k) + s + 1 );
         *( *(uijx + k) + s ) = *( *(uijx + k) + s + 1 );
         *( *(uijy + k) + s ) = *( *(uijy + k) + s + 1 );
         *( *(uijNSQ + k) + s ) = *( *(uijNSQ + k) + s + 1 );
//         *( *(aijASQ + k) + s ) = *( *(aijASQ + k) + s + 1 );
         *( *(Ouijx + k) + s ) = *( *(Ouijx + k) + s + 1 );
         *( *(Ouijy + k) + s ) = *( *(Ouijy + k) + s + 1 );
             
    }
         
    fn[k] = fn[k] - 1;
         
    for (int r = 0; r < fn[bk]; r++) {      // Find out bj.        
        
        if ( *(*( nnl + bk ) + r ) == k ) {        
             bj = r;
             break;             
        }
        
    }
    
    if (bj == -1) {cout << "bj is wrong!" << endl;} 
             
    for (int s = bj; s < (fn[bk]-1); s++) {
             
         *( *(nnl + bk) + s ) = *( *(nnl + bk) + s + 1 );
         *( *(uijx + bk) + s ) = *( *(uijx + bk) + s + 1 );
         *( *(uijy + bk) + s ) = *( *(uijy + bk) + s + 1 );
         *( *(uijNSQ + bk) + s ) = *( *(uijNSQ + bk) + s + 1 );
//         *( *(aijASQ + bk) + s ) = *( *(aijASQ + bk) + s + 1 );
         *( *(Ouijx + k) + s ) = *( *(Ouijx + k) + s + 1 );
         *( *(Ouijy + k) + s ) = *( *(Ouijy + k) + s + 1 );
             
    }
             
    fn[bk] = fn[bk] - 1;

}

void uijComponent (int** ListH, int* LengthH, double* ArrayH, int NP, double** TBWritten) {     // (nnl, &fn[0], &ux[0], NP, uijx)

    for (int k = 0; k < NP; k++) {
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            *( TBWritten[k] + j ) = *( ArrayH + *( *(ListH + k) + j ) ) - *( ArrayH + k );
            
        }
        
    }

}

void OuijCC (double** uijCO, int* LengthH, int NP, double** TBWritten, double Stch) {     // Arg: (uijx, &fn[0], NP, Ouijx, StchX)

    for (int k = 0; k < NP; k++) {
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            *( TBWritten[k] + j ) = ( *( *(uijCO + k) + j ) )/Stch ;
            
        }
        
    }

}
/*
void aijdotASQ (int** NnlH, int* LengthH, double** OuijCHX, double** OuijCHY, int NP, double Ax, double Ay, double** TBW, double stchX, double stchY) {    // Arg: (nnl, &fn[0], uijx, uijy, NP, Ax, Ay, aijASQ, StchX, StchY)
     
     for (int k = 0; k < NP; k++) {
         
         for (int j = 0; j < *(LengthH + k); j++) { 
                      
             *( TBW[k] + j) = pow( *(*(OuijCHX + k) + j)*Ax/stchX + *(*(OuijCHY + k) + j)*Ay/stchY , 2);   
                   
         }
         
     }
     
}
*/
void uijNormSQ (int** nnl, double** uijCH1, double** uijCH2, double** Ouijx, double** Ouijy, int* LengthH, int NP, double** TBWritten, double* K, int* broken, int n, int time_no_snap) {              // Arg: (nnl, uijx, uijy, aijASQ, &fn[0], NP, uijNSQ, K)

    double bff = 0.0;
    
    for (int k = 0; k < NP; k++) {
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
             
             bff = *( *(uijCH1 + k) + j ) * *( *(uijCH1 + k) + j )   
                      + *( *(uijCH2 + k) + j ) * *( *(uijCH2 + k) + j );
                                  
             if ( n >= time_no_snap ) { *( TBWritten[k] + j ) = bff; }
             
             else if( bff <= (Lf0+(g/K[k]))*(Lf0+(g/K[k]))*a*a ) {             
                 *( TBWritten[k] + j ) = bff;           
             }   
             
             else {
                  UpdateList(nnl, uijCH1, uijCH2, Ouijx, Ouijy, TBWritten, LengthH, NP, k, j);
                  broken[k] = broken[k] + 1;     
             }
             
        }
        
    }
    
}

void FuComponent (double** uijCH, int* LengthH, double a, int NP, double* TBWritten, double** uijnormsq) {             // Arg: (uijx, &fn[0], a, NP, xFu, uijNormSQ)
        
       for (int k = 0; k < NP; k++) {
        
        TBWritten[k] = 0.0;
        
        double ratio = *(LengthH + k)/6.0;
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            TBWritten[k] = TBWritten[k] + *(*(uijCH + k) + j);
            
        }
        
        TBWritten[k] = (-1)*TBWritten[k]/(3*ratio*a*a);
        
    }
    
}
/*
void KuComponent (double** uijCH1, int* LengthH, double a, int NP, double** uijCH2, double* TBWritten) {       
                                                                                                       // Arg: (uijx, &fn[0], a, NP, uijy, xKu)
    
    for (int k = 0; k < NP; k++) {
        
        TBWritten[k] = 0.0;
        
        double ratio = (*(LengthH + k)/6.0) * (*(LengthH + k)/6.0);
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            for (int m = 0; m < *(LengthH + k); m++) {
            
                if ( j != m ) {
                          
                    TBWritten[k] = TBWritten[k] + ( ( *(*(uijCH1 + k) + j) )*( *(*(uijCH2 + k) + m) ) - ( *(*(uijCH1 + k) + m) )*( *(*(uijCH2 + k) + j) ) ) *  ( *(*(uijCH2 + k) + j)  - *(*(uijCH2 + k) + m) );
                
                }
            }
        
        }
        
        TBWritten[k] = TBWritten[k]/(9*ratio*a*a*a*a);
        
    }

}  
*/
/*
void MuComponent (double** uijCH1, int* LengthH, double a, int NP, double* AComp1, double* AComp2, double** uijCH2, double* TBW) {     // Arg: (uijx, &fn[0], a, NP, &Ax[0], &Ay[0], uijy, xMu)

     for (int k = 0; k < NP; k++) {
         
         TBW[k] = 0.0;
         
         double ratio = *(LengthH + k)/6.0;
         
         for (int j=0; j < *(LengthH + k); j++) {
             
             TBW[k] = TBW[k] + AComp1[k]*( *(uijCH1[k] + j)*AComp1[k] + *(uijCH2[k] + j)*AComp2[k] );
             
         }
              
         TBW[k] = TBW[k]/(-3*ratio*a*a);
         
     }
      
}
*/

void Fsum (double** uijnormsq, int* LengthH, double a, int NP, double* TBWritten, int* broken, double* K) {                // Arg: (uijNSQ, &fn[0], a, NP, Fsum)
    
    for (int k = 0; k < NP; k++) {
        
        TBWritten[k] = 0.0;
        
        double ratio = *(LengthH + k)/6.0;
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            TBWritten[k] = TBWritten[k] + *(*(uijnormsq + k) + j) - (a*a);
            
        }
        
        TBWritten[k] = TBWritten[k] + ( (Lf0+(g/K[k]))*(Lf0+(g/K[k]))*a*a - (a*a) )*broken[k];
        
        TBWritten[k] = TBWritten[k]/(6*ratio*a*a);
        
    }
  
}

void Ksum(double** uijCH1, int* LengthH, double a, int NP, double** uijCH2, double* TBWritten) {       
                                                                                         // Arg: (uijx, uijNSQ, &fn[0], a, NP, uijy, Ksum)
    
    for (int k = 0; k < NP; k++) {
        
        TBWritten[k] = 0.0;
        
        double ratio = (*(LengthH + k)/6.0) * (*(LengthH + k)/6.0);
        
        for (int j = 0; j < *(LengthH + k); j++) {
            
            for (int m = 0; m < *(LengthH + k); m++) {
            
                if ( j != m ) {
                                 
                    TBWritten[k] = TBWritten[k] + ( ( *(*(uijCH1 + k) + j) )*( *(*(uijCH2 + k) + m) ) - ( *(*(uijCH1 + k) + m) )*( *(*(uijCH2 + k) + j) ) ) 
                                                  * ( ( *(*(uijCH1 + k) + j) )*( *(*(uijCH2 + k) + m) ) - ( *(*(uijCH1 + k) + m) )*( *(*(uijCH2 + k) + j) ) ) ; 
                    
                }
                
                
            }
        
        }
        
        TBWritten[k] = TBWritten[k]/(18*ratio*a*a*a*a);
        
    }
    
}  

void Msum (double** uijCHX, double** uijCHY, double** Ouijx, double** Ouijy, int* LengthH, int NP, double a, double Ax, double Ay, double* TBW) {       // Arg: (uijx, uijy, Ouijx, Ouijy, &fn[0], NP, a, &Ax[0], &Ay[0], M)

     for (int k = 0; k < NP; k++) {
         
         TBW[k] = 0.0;
         
         double ratio = *(LengthH + k)/6.0;
         
         for (int j = 0; j < *(LengthH + k); j++) {
             
             TBW[k] = TBW[k] +  pow( *(uijCHX[k] + j)*Ax + *(uijCHY[k] + j)*Ay , 2) - pow( *(Ouijx[k] + j)*Ax + *(Ouijy[k] + j)*Ay , 2);
             
         }
         
         TBW[k] = TBW[k]/(6*ratio*a*a);
         
     }    
     
}


double GoUValue_CL (double F, double K, double M, double FuComp, double KuComp, double MuComp, double m, double A, double B, double valley) {           //Just for calculation

    double UI1 = m*A; 
    double UI2 = 0;
    double UI4 = 0;
   
    return UI1 * FuComp + UI2 * (0.25 * KuComp - 0.5 * FuComp) + UI4 * MuComp;
    
}
/*
double GoUValue_RX (double F, double K, double M, double FuComp, double KuComp, double MuComp, double m, double A, double B, double valley) {           //Just for calculation

    double UI1 = 2*0.18*m*A*( 0.00125*(2*F+3)*(2*F+3)*(2*F+3) - 0.07*(2*F+3)*(2*F+3) + (2*F+3) ); 
    double UI2 = 0;
    double UI4 = 0;
   
    return UI1 * FuComp + UI2 * (0.25 * KuComp - 0.5 * FuComp) + UI4 * MuComp;
    
} 
*/
void arrtr (double* original, double* written, int NP){          
                                               //[Array Transcript] Arg: eg. (&uxnew[0], &uxold[0], NP) -> uxold will be written by uxnew.
    
    for (int k = 0; k < NP; k++) {
        
        *(written + k) = * (original + k);
        
    }
    
}

void updA (double** uijx, double** uijy, double** uijNSQ, int* LengthH, int* AFF, double a, double CRIT, double* Ax, double* Ay) {    // Arg: (k, uijx, uijy, uijNSQ, &AFF[0], a, CRIT, &fn[0], &Ax[0], &Ay[0])

          for (int k = 0; k < NP; k++) {
              
              if ( AFF[k] != 1 ) { continue; }
              
              else {              
                       
                  double xc = 0;
                  double yc = 0; 
              
                  for (int j = 0; j < *(LengthH + k); j++ ){
                      
                      if ( (*(*(uijx + k) + j))*(*(*(uijy + k) + j)) > 0 ) {          // The "equal to" sign is protected by the AFF=2 for original setting.
                           
                           if ( *(*(uijNSQ + k) + j) >= CRIT*CRIT*a*a ) {
                               xc = xc + abs(*(*(uijx + k) + j));
                               yc = yc + abs(*(*(uijy + k) + j));
//                               if ((n==19) && (k==tork)) {cout << "{" << abs(*(*(uijx + k) + j)) << ", " << abs(*(*(uijy + k) + j)) << ")" << endl; }
                           }
                      }
                      
                      if ( (*(*(uijx + k) + j))*(*(*(uijy + k) + j)) < 0 ) {
                           
                           if ( (*(*(uijx + k) + j)) < 0 ) {               // Not yet considered the equal sign's possible problem
                                
                              if ( *(*(uijNSQ + k) + j) >= CRIT*CRIT*a*a ) {
                                   xc = xc + (*(*(uijx + k) + j));
                                   yc = yc + (*(*(uijy + k) + j));
//                                   if ((n==19) && (k==tork)) {cout << "{" << (*(*(uijx + k) + j)) << ", " <<  (*(*(uijy + k) + j)) << ")" << endl; }
                              }
                              
                           }
                           
                           if ( (*(*(uijx + k) + j)) > 0 ) {
                               
                               if ( *(*(uijNSQ + k) + j) >= CRIT*CRIT*a*a ) {
                                   xc = xc + (-1)*(*(*(uijx + k) + j));
                                   yc = yc + (-1)*(*(*(uijy + k) + j));
//                                   if ((n==19) && (k==tork)) {cout << "{" << (-1)*(*(*(uijx + k) + j)) << ", " << (-1)*(*(*(uijy + k) + j)) << ")" << endl; }
                              }
                                                             
                           }
                           
                      }
                                   
                  }
                  
                  if ( (xc == 0) && (yc == 0) ) { continue; }                
                  *(Ax + k) = xc/sqrt(xc*xc + yc*yc); 
                  *(Ay + k) = yc/sqrt(xc*xc + yc*yc); 
//                  if ((n==19) && (k==tork)) {cout << "{" << xc << ", " << yc << ") " << sqrt(xc*xc + yc*yc) << endl; }
                  
              }
              
          }
          
}

//
int main () {
    
    double ux[NP];
    double uy[NP];
    const double InX = 1.0*a;                   // Increment of X direction
    const double InY = 2*0.8660254*a;               // Increment of Y direction
    int boxi[NP];
    int boxj[NP];
    int lbl = 0;                    // List of particle; particle label; initialization
    int wx = static_cast<int>( (0.5*X*a+0.5*X*a+0.25*a)/InX); 
    int wy = static_cast<int>( abs(-0.8660254*Y*a-0.8660254*Y*a-0.5*0.8660254*a)/InY); 
    int nob = (wx+1)*(wy+1);                         // number of box
    int BLL = 2;                              // Box List Length
    int eta[NP];
    int kappa[NP];
//    double Ax[NP]; 
//    double Ay[NP];
    int AFF[NP];
    double valley[NP];

    for (int k = 0; k < NP; k++) {
        
//        Ax[k] = 0;
//        Ay[k] = 1;
        AFF[k] = 2;
        valley[k] = onset;    
        
    }
    
    for (int h = -1*X; h <= X; h++){
   
        for (int k = 0; k <= Y; k++){
            
            if ( abs(0.5*a*h + 0.5*a*k) > X*0.5*a ) {
                break;
            }
                
            else if (abs(0.8660254*a*h - 0.8660254*a*k) > X*0.5*a) {
                continue;
            }
            
            else
            ux[lbl] = 0.5*a*h + 0.5*a*k;           // u = h*a + k*b
            uy[lbl] = 0.8660254*a*h - 0.8660254*a*k;       // a = (1/2, sqrt(3)/2), b = (1/2, -sqrt(3)/2)
            boxi[lbl] = static_cast<int>((ux[lbl]+0.5*X*a+0.25*a)/InX);
            boxj[lbl] = static_cast<int>(abs(uy[lbl]-0.8660254*Y*a-0.5*0.8660254*a)/InY);
            eta[lbl] = h;
            kappa[lbl] = k;
            lbl++;
            
        }    
     
        for (int k = -1; k >= -1*Y; k--){
            
            if ( abs(0.5*a*h + 0.5*a*k) > X*0.5*a ) {
                break;
            }
                
            else if (abs(0.8660254*a*h - 0.8660254*a*k) > X*0.5*a) {
                continue;
            }
            
            else
            ux[lbl] = 0.5*a*h + 0.5*a*k;           // u = h*a + k*b
            uy[lbl] = 0.8660254*a*h - 0.8660254*a*k;       // a = (1/2, sqrt(3)/2), b = (1/2, -sqrt(3)/2)
            boxi[lbl] = static_cast<int>((ux[lbl]+0.5*X*a+0.25*a)/InX);
            boxj[lbl] = static_cast<int>(abs(uy[lbl]-0.8660254*Y*a-0.5*0.8660254*a)/InY);
            eta[lbl] = h;
            kappa[lbl] = k;
            lbl++;
            
        }        
          
    }    

    int **boxlist = (int**) malloc (nob*sizeof(int*)) ;
    int *boxptr = (int*) malloc (nob*BLL*sizeof(int)) ;
    int *aux = (int*) calloc (nob, sizeof(int));  
    
    for (int i = 0; i < nob; i++) {
        
        boxlist[i] = boxptr;
        boxptr += BLL;
    
    }
    
    for (int k = 0; k < NP; k++) {
        
        int p = (wx+1) * boxj[k] + boxi[k];
        *(boxlist[p] + aux[p]) = k;
        aux[p]++;
        
    }
    
    int **nnl = (int**) malloc (NP*sizeof(int*));          // Nearest Neighbor List (nnl is the first pointer of pointers)
    int *fn = (int*) calloc (NP, sizeof(int));             // First Neighbor (array)
    int *broken = (int*) calloc (NP, sizeof(int));
    int adj[9] = {0, -1, 1, -1*(wx+1), -1-(wx+1), 1-(wx+1), (wx+1), -1+(wx+1), 1+(wx+1)};

    
    for (int k = 0; k < NP; k++) {                   // In order to find fn array
        
        int bxl = (wx+1)*boxj[k] + boxi[k];        // the particle's box number
        int shxtmp;                                // shxtmp means the box being focused on
        double lengthSQ;    
           
        for (int j = 0; j < 9; j++) {                       // Is there any more general function to insert 5, the size of the char array?
                    
            shxtmp = bxl + adj[j] ;      // problematic! (Solved)
            if (shxtmp < 0 || shxtmp >= nob) {              // nob is the number of box
                continue;
            }
               
            for (int i = 0; i < aux[shxtmp]; i++) {             // The upper bound should be aux[shxtmp], not BLL.
                
                lengthSQ = ENSQ(ux[k], uy[k], ux[*(boxlist[shxtmp]+i)], uy[*(boxlist[shxtmp]+i)]);    
                if (0.99*a*a < lengthSQ && lengthSQ < 1.01*a*a) {
                    fn[k]++;                                                                // fn is not counting the number!
               }
                          
            }
            
        }        
        
    }
    
    for (int k = 0; k < NP; k++) {
        
        int* nnptr = (int*) malloc(fn[k]*sizeof(int));
        nnl[k] = nnptr;
        
    }
    
    for (int k = 0; k < NP; k++) {                 // Really assign the particle numbers into nnl
        

        int bxl = (wx+1)*boxj[k] + boxi[k];        // the particle's box number
        int shxtmp;                                // shxtmp means the box being focused on
        double lengthSQ;
        int count = 0;    
           
        for (int j = 0; j < 9; j++) {                       // Is there any more general function to insert 5, the size of the char array?
                    
            shxtmp = bxl + adj[j] ;      // problematic!
            if (shxtmp < 0 || shxtmp >= nob) {              // nob is the number of box
                continue;
            }
               
            for (int i = 0; i < aux[shxtmp]; i++) {             // The upper bound should be aux[shxtmp], not BLL.
                
                lengthSQ = ENSQ(ux[k], uy[k], ux[*(boxlist[shxtmp]+i)], uy[*(boxlist[shxtmp]+i)]);    
                if (0.99*a*a < lengthSQ && lengthSQ < 1.01*a*a) {
                    *(nnl[k]+count) = *(boxlist[shxtmp]+i);
                    count++;  
                    if (count > fn[k]) {
                        cout << "WRONG!" << endl;    
                    }                                               
                }
                          
                
                
            }
            
        }        
        
    } 
 
    for ( int k = 0; k < NP; k++ ) {
        
        uy[k] = stchY*uy[k];
        ux[k] = stchX*ux[k];
        
    } 
 
    int count1 = 0;  // long
    int count2 = 0;  // short
    int *lg = (int*) calloc (SLL, sizeof(int));
    int *sht = (int*) calloc (SLL-1, sizeof(int));


    for (int k = 0; k < NP; k++) {
         
         if ( eta[k]-kappa[k] == 0 && ux[k] <= (ER + 0.1)*stchX*a && ux[k] >= (EL - 0.1)*stchX*a ) {
              
              lg[count1] = k;
              count1++;
              
         }
        
    }

//    cout << endl;

    if (count1 != SLL) {cout << "Wrong1!" << endl; }

    for (int k = 0; k < NP; k++) {
         
         if ( eta[k]-kappa[k] == 1 && (ux[k] < (ER - 0.1)*stchX*a) && (ux[k] > (EL + 0.1)*stchX*a ) ) {
              
              sht[count2] = k;
              count2++;
                           
         }
        
    }

//    cout << endl;

    if (count2 != SLL-1) {cout << "Wrong2!" << endl; }

    for (int mm = 0; mm < SLL; mm++) {
        
        for (int n = 0; n < SLL-1; n++) {
            
            for (int r = 0; r < fn[lg[mm]]; r++) {    // loop for matching the upper rim particle (sht[n]).
                
                if (  *(nnl[lg[mm]] + r ) == sht[n] ) {
                     
                     for ( int s = r; s < (fn[lg[mm]]-1); s++) {
                         
                         *( nnl[lg[mm]] + s ) = *( nnl[lg[mm]] + s + 1 );
                         
                     } 
                     
                     fn[lg[mm]] = fn[lg[mm]] - 1;
                     
                }
                
            }
                
        }
    
    }
    
    for (int mm = 0; mm < SLL-1; mm++) {
    
        for (int n = 0; n < SLL; n++) {    
            
            for (int r = 0; r < fn[sht[mm]]; r++) {    // loop for matching the upper rim particle (sht[n]).
                
                if ( *(nnl[sht[mm]] + r ) == lg[n] ) {
                     
                     for ( int s = r; s < (fn[sht[mm]]-1); s++) {
                         
                         *( nnl[sht[mm]] + s ) = *( nnl[sht[mm]] + s + 1 );
                         
                     } 
                     
                     fn[sht[mm]] = fn[sht[mm]] - 1;
                     
                }
                
            }
                
        }
    
    }
/* 
    for (int j = 0; j < SLL-1; j++) {
        
        cout << sht[j] << '(' << fn[sht[j]] << ')' << ((fn[sht[j]] == 5) ? "* " : " ");
    
    }
    
    cout << endl;
     
    for (int j = 0; j < SLL; j++) {
        
        cout << lg[j] << '(' << fn[lg[j]] << ')' << ((fn[lg[j]] == 5) ? "* " : " ");
    
    }
    
    cout << endl << endl;
*/
    double **uijx = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's uijx   

    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        uijx[k] = nnptr;
        
    }
    
    double **uijy = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's uijy
    
    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        uijy[k] = nnptr;
        
    }
    
    double **uijNSQ = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's uij norm squared
    
    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        uijNSQ[k] = nnptr;
        
    }
/*
    double **aijASQ = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's uij norm squared
    
    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        aijASQ[k] = nnptr;
        
    }
*/
    double **Ouijx = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's Ouijx   

    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        Ouijx[k] = nnptr;
        
    }

    double **Ouijy = (double**) malloc (NP*sizeof(double*));          // Nearest Neighbor List's Ouijy
    
    for (int k = 0; k < NP; k++) {
        
        double* nnptr = (double*) malloc(fn[k]*sizeof(double));
        Ouijy[k] = nnptr;
        
    }
    
    double* F = (double*) malloc (NP* sizeof(double));
    double* K = (double*) malloc (NP* sizeof(double));
    
    double* Mx = (double*) malloc (NP* sizeof(double));
    double* My = (double*) malloc (NP* sizeof(double));
    
    double* xFu = (double*) malloc (NP* sizeof(double));
    double* yFu = (double*) malloc (NP* sizeof(double));
    double* xKu = (double*) malloc (NP* sizeof(double));
    double* yKu = (double*) malloc (NP* sizeof(double));
    double* xMu = (double*) malloc (NP* sizeof(double));
    double* yMu = (double*) malloc (NP* sizeof(double));   
    
    for (int k = 0; k < NP; k++) {
     
        K[k] = K0;
        
    }  
    
    double* GoUx = (double*) malloc (NP* sizeof(double));
    double* GoUy = (double*) malloc (NP* sizeof(double));
    double* PastGoUx = (double*) malloc (NP* sizeof(double));
    double* PastGoUy = (double*) malloc (NP* sizeof(double)); 
    double* uddx = (double*) malloc (NP* sizeof(double));
    double* uddy = (double*) malloc (NP* sizeof(double));
    double* uxnew = (double*) malloc (NP* sizeof(double));
    double* uxold = (double*) malloc (NP* sizeof(double));
    double* uynew = (double*) malloc (NP* sizeof(double));
    double* uyold = (double*) malloc (NP* sizeof(double));    

    ////////////////////////////0723 relax/////////////////////////

    int* relax = (int*) malloc (NP* sizeof(int));
    double* ETA = (double*) malloc (NP* sizeof(double));

    for (int k = 0; k < NP; k++) {

      relax[k] = 1;

    }

    ////////////////////////////0723 relax//////////////////////////    

    uijComponent(nnl, &fn[0], &ux[0], NP, uijx);
    uijComponent(nnl, &fn[0], &uy[0], NP, uijy);    
//    aijdotASQ(nnl, &fn[0], uijx, uijy, NP, Ax, Ay, aijASQ, stchX, stchY);
    OuijCC(uijx, &fn[0], NP, Ouijx, stchX);         // Calculate for Ouijx and Ouijy    
    OuijCC(uijy, &fn[0], NP, Ouijy, stchY);
    uijNormSQ(nnl, uijx, uijy, Ouijx, Ouijy, &fn[0], NP, uijNSQ, K, broken, 0, time_no_snap);


    FuComponent(uijx, &fn[0], a, NP, xFu, uijNSQ);
    FuComponent(uijy, &fn[0], a, NP, yFu, uijNSQ);
    
//    KuComponent(uijx, &fn[0], a, NP, uijy, xKu);
//    KuComponent(uijy, &fn[0], a, NP, uijx, yKu);
    
//    MuComponent(uijx, &fn[0], a, NP, &Ax[0], &Ay[0], uijy, xMu);
//    MuComponent(uijy, &fn[0], a, NP, &Ay[0], &Ax[0], uijx, yMu);    // Note that the order of y and x are special!
    
    Fsum(uijNSQ, &fn[0], a, NP, F, broken, K);
    Ksum(uijx, &fn[0], a, NP, uijy, K);


    for (int k = 0; k < NP; k++) {             
                 
                GoUx[k] = GoUValue_CL(F[k], K[k], My[k], xFu[k], xKu[k], xMu[k], m, A, B, valley[k]);
      
                GoUy[k] = GoUValue_CL(F[k], K[k], My[k], yFu[k], yKu[k], yMu[k], m, A, B, valley[k]);
     
    }
/*    
        cout << "T = -1" << endl  
             << "k = " << tork << endl
             << "fn = " << fn[tork] << endl << endl
             << "AFF = " << AFF[tork] << endl << endl
             << "ux = " << ux[tork] << endl
             << "uy = " << uy[tork] << endl << endl
             << "Ax = " << Ax[tork] << endl
             << "Ay = " << Ay[tork] << endl << endl
             << "eta = " << eta[tork] << endl
             << "kappa = " << kappa[tork] << endl << endl
             << "a = " << a << endl << endl
             << endl << "--------" << endl << endl;
                                       
        shu = fn[tork];
             
        cout << "Neighbor vectors (uij's) : " << endl;
             
        for (int j = 0; j < fn[tork]; j++) {
            
            cout<< "(" << (*(*(uijx + tork)+j))/a << ", " << (*(*(uijy + tork)+j))/a << ")" << endl;     
            
        }
    
        cout << "--------" << endl;
*/    
    arrtr(&ux[0], &uxold[0], NP);
    arrtr(&uy[0], &uyold[0], NP);
    


    
    for (int k = 0; k < NP; k++) {
        
        if ( (eta[k]+kappa[k] == -X) || (eta[k]+kappa[k] == -X+1) ) { uddx[k] = 0; uddy[k] = 0; }
        else if ( (eta[k]+kappa[k] == X) || (eta[k]+kappa[k] == X-1) ) { uddx[k] = 0; uddy[k] = 0; } 
        else if ( eta[k]-kappa[k] == (2*floor(0.5*X/1.7320508)+1) ) { uddx[k] = 0; uddy[k] = 0; }
        else if ( eta[k]-kappa[k] == -(2*floor(0.5*X/1.7320508)+1) ) { uddx[k] = 0; uddy[k] = 0; } 
        
        else {
        uddx[k] = (1/m)*(-GoUx[k]);                        // Initiate the acceleration in x direction
        uddy[k] = (1/m)*(-GoUy[k]);                        // Initiate the acceleration in y direction
        }
        
        ux[k] = ux[k] + 0.5*uddx[k]*dt*dt;                               
        uy[k] = uy[k] + 0.5*uddy[k]*dt*dt;
        
               
    }
/*   
        cout << "T = 0" << endl  
             << "k = " << tork << endl
             << "fn = " << fn[tork] << endl << endl
             << "AFF = " << AFF[tork] << endl << endl
             << "ux = " << ux[tork] << endl
             << "uy = " << uy[tork] << endl << endl
             << "Ax = " << Ax[tork] << endl
             << "Ay = " << Ay[tork] << endl << endl
             << "F = " << F[tork] << endl
             << "K = " << K[tork] << endl
             << "M = " << M[tork] << endl << endl
             
             << "xFu = " << xFu[tork] << endl 
             << "yFu = " << yFu[tork] << endl << endl
             
             << "xKu = " << xKu[tork] << endl
             << "yKu = " << yKu[tork] << endl << endl
             
             << "xMu = " << xMu[tork] << endl
             << "yMu = " << yMu[tork] << endl << endl
             
             << "-GoUx = " << -1*GoUx[tork] << endl
             << "-GoUy = " << -1*GoUy[tork] << endl
             << "-xDissi = " << -beta*(GoUx[tork]-PastGoUx[tork])/dt << endl
             << "-yDissi = " << -beta*(GoUy[tork]-PastGoUy[tork])/dt << endl 
             << endl << "--------" << endl << endl;
*/
   
    arrtr(&ux[0], &uxnew[0], NP);
    arrtr(&uy[0], &uynew[0], NP);
    arrtr(&GoUx[0], &PastGoUx[0], NP);
    arrtr(&GoUy[0], &PastGoUy[0], NP);

// Iteration of time

    for (int n = 1; n <= T; n++) {
        

// Update based on the last ux[k] and uy[k].        
        uijComponent(nnl, &fn[0], &ux[0], NP, uijx);                     // The n=1 is the first time to use the changed uij's to calculate the numerical functions. 
        uijComponent(nnl, &fn[0], &uy[0], NP, uijy);    
        uijNormSQ(nnl, uijx, uijy, Ouijx, Ouijy, &fn[0], NP, uijNSQ, K, broken, n, time_no_snap);       // No aijdotASQ because it shouldn't change with time.      
    
//        updA(uijx, uijy, uijNSQ, &fn[0], &AFF[0], a, CRIT, &Ax[0], &Ay[0]);
    
        FuComponent(uijx, &fn[0], a, NP, xFu, uijNSQ);
        FuComponent(uijy, &fn[0], a, NP, yFu, uijNSQ);
        
//        KuComponent(uijx, &fn[0], a, NP, uijy, xKu);
//        KuComponent(uijy, &fn[0], a, NP, uijx, yKu);
        
//        MuComponent(uijx, &fn[0], a, NP, &Ax[0], &Ay[0], uijy, xMu);
//        MuComponent(uijy, &fn[0], a, NP, &Ay[0], &Ax[0], uijx, yMu);
        
        Fsum(uijNSQ, &fn[0], a, NP, F, broken, K);
        Ksum(uijx, &fn[0], a, NP, uijy, K);

        
        
        for (int k = 0; k < NP; k++) {
                 
                GoUx[k] = GoUValue_CL(F[k], K[k], My[k], xFu[k], xKu[k], xMu[k], m, A, B, valley[k]);
      
                GoUy[k] = GoUValue_CL(F[k], K[k], My[k], yFu[k], yKu[k], yMu[k], m, A, B, valley[k]);

        }
        
        
        for (int k = 0; k < NP; k++) {
        
            if ( (eta[k]+kappa[k] == -X) || (eta[k]+kappa[k] == -X+1) ) { uddx[k] = 0; uddy[k] = 0; }
            else if ( (eta[k]+kappa[k] == X) || (eta[k]+kappa[k] == X-1) ) { uddx[k] = 0; uddy[k] = 0; } 
            else if ( eta[k]-kappa[k] == (2*floor(0.5*X/1.7320508)+1) ) { uddx[k] = 0; uddy[k] = 0; }
            else if ( eta[k]-kappa[k] == -(2*floor(0.5*X/1.7320508)+1) ) { uddx[k] = 0; uddy[k] = 0; } 
            
            else {
            uddx[k] = (1/m)*(-GoUx[k]-beta*(GoUx[k]-PastGoUx[k])/dt);                        // Initiate the acceleration in x direction
            uddy[k] = (1/m)*(-GoUy[k]-beta*(GoUy[k]-PastGoUy[k])/dt);                        // Initiate the acceleration in y direction
            }
            ux[k] = 2*uxnew[k] - uxold[k] + uddx[k]*dt*dt;
            uy[k] = 2*uynew[k] - uyold[k] + uddy[k]*dt*dt;
            
//            if ( abs(uy[k]) > 0.5*sqrt(3)*a*X ) { cout << "Exceeded!!" << endl; }
            
        }
/*
        if ( (n%100 == 0) || (fn[tork]!=shu) ) {
                  
            cout << "T = " << n << endl  
                 << "k = " << tork << endl
                 << "fn = " << fn[tork] << endl << endl
                 << "AFF = " << AFF[tork] << endl << endl
                 << "ux = " << ux[tork] << endl
                 << "uy = " << uy[tork] << endl << endl
                 << "Ax = " << Ax[tork] << endl
                 << "Ay = " << Ay[tork] << endl << endl
                 << "F = " << F[tork] << endl
                 << "K = " << K[tork] << endl
                 << "M = " << M[tork] << endl << endl
                 
                 << "xFu = " << xFu[tork] << endl 
                 << "yFu = " << yFu[tork] << endl << endl
                 
                 << "xKu = " << xKu[tork] << endl
                 << "yKu = " << yKu[tork] << endl << endl
                 
                 << "xMu = " << xMu[tork] << endl
                 << "yMu = " << yMu[tork] << endl << endl
                 
                 << "Ax = " << Ax[tork] << endl
                 << "Ay = " << Ay[tork] << endl << endl
                 
                 << "-GoUx = " << -1*GoUx[tork] << endl
                 << "-GoUy = " << -1*GoUy[tork] << endl
                 << "-xDissi = " << -beta*(GoUx[tork]-PastGoUx[tork])/dt << endl
                 << "-yDissi = " << -beta*(GoUy[tork]-PastGoUy[tork])/dt << endl << endl; 
                 
                 for ( int j = 0; j < fn[tork]; j++) {
      
                     cout << *(nnl[tork]+j) << " (" << sqrt( *(uijNSQ[tork]+j) )/a << ')' 
                          << " [" << *(uijx[tork]+j) /a << ", " << *(uijy[tork]+j) /a << ']' << endl;
                          
                 }  
                 
                 cout << endl << "--------" << endl << endl;
                 
        }
        
        if (fn[tork]!=shu) { shu = fn[tork]; }                           
 
*/           
        arrtr(&uxnew[0], &uxold[0], NP);
        arrtr(&uynew[0], &uyold[0], NP);
        arrtr(&ux[0], &uxnew[0], NP);
        arrtr(&uy[0], &uynew[0], NP);
        arrtr(&GoUx[0], &PastGoUx[0], NP);
        arrtr(&GoUy[0], &PastGoUy[0], NP);
      
        if ( (n%50000==0) || (n==T) || (n==T-1) ||  (n==3000) ) {

                 ofstream fout ( (  string("NH_noM_dt4e-7")
                                  + string("_T")
                                  + int_to_string(n)
                                  + string(".txt")       ).c_str() 
                                 );
                 
                 for ( int k = 0; k < NP; k++ ) {
                     
                     fout << ux[k] << ' ' << uy[k] << endl;
                     
                 }
                             
             } // end of if ( (n%1000 == 0) || (n == 1) ) 
            



    } // end of n loop
    
        Msum(uijx, uijy, Ouijx, Ouijy, &fn[0], NP, a, 1, 0, Mx);
        Msum(uijx, uijy, Ouijx, Ouijy, &fn[0], NP, a, 0, 1, My);
        
        ofstream  Mx_out   ( (  string("NH_noM_dt4e-7")
                                  + string("_T")
                                  + int_to_string(T)
                                  + string("Mx.txt")       ).c_str() 
                                 );
                                 
        ofstream  My_out   ( (  string("NH_noM_dt4e-7")
                                  + string("_T")
                                  + int_to_string(T)
                                  + string("My.txt")       ).c_str() 
                                 );
        
        for ( int k = 0; k < NP; k++ ) {
                     
                Mx_out << ux[k] << ' ' << uy[k] << ' ' << sqrt(2*Mx[k]+1) << endl;
                My_out << ux[k] << ' ' << uy[k] << ' ' << sqrt(2*My[k]+1) << endl; 
                     
        } // end of for loop over k
    



    ofstream nnl_output( ( string("NH_dt2e-8") + string("nnl_T") + int_to_string(T) + string(".txt") ).c_str() );
    ofstream fn_output( ( string("NH_dt2e-8") + string("fn_T") + int_to_string(T) + string(".txt") ).c_str() );
    
    for ( int k = 0; k < NP; k++ ) {
        
        fn_output << fn[k] << endl;
        
    }
    
    for ( int k = 0; k < NP; k++ ) {
        
        for ( int j = 0; j < fn[k]; j++ ) {
  
            nnl_output << *(*(nnl+k)+j) << " ";
            
        }
        
        nnl_output << '\n';
        
    }

  

//system("pause");
    
}
