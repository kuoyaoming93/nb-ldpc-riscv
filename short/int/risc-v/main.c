#include <stdint.h>
#include "printf.h"
#include "runtime.h"


#include "constantes_CONF.h"
#include "gf_tables_GF16.h"
#include "CNU_tables_GF16.h"
#include "Hmat_N32_M16_GF16.h"

//#define VECTOR_EXT
//#define FOR_BACK

int16_t Ln_aux[16][32] = {
    {6,2,7,0,5,0,3,2,11,4,0,1,6,5,3,4,6,7,9,1,5,12,2,1,12,3,3,2,2,4,4,5},
    {0,0,9,0,4,2,6,3,7,7,3,5,2,2,6,6,5,6,6,0,3,9,3,3,8,0,4,6,6,2,8,1},
    {9,1,4,3,3,0,0,4,6,2,2,3,5,8,4,0,10,4,5,4,4,9,1,0,8,4,5,5,1,2,3,4}, 
    {8,2,5,0,9,3,3,2,9,2,1,2,4,8,2,7,3,6,6,5,7,8,2,4,13,3,1,0,6,7,4,9}, 
    {6,6,5,0,3,3,6,0,10,4,0,0,8,3,0,7,3,4,10,1,3,10,1,2,7,8,1,6,0,6,2,4}, 
    {3,0,6,3,2,1,3,5,2,5,4,7,1,5,7,1,9,3,2,3,2,6,1,1,4,0,6,9,5,0,6,1}, 
    {11,2,1,3,7,2,0,4,4,0,3,4,3,11,3,2,8,3,2,7,5,5,0,2,9,3,3,2,5,4,2,8}, 
    {8,7,3,1,7,6,6,0,8,2,1,1,6,6,0,10,1,3,7,4,5,6,1,5,8,8,0,4,4,9,1,8}, 
    {2,4,4,3,0,4,6,3,1,5,4,5,3,3,4,4,7,0,3,2,0,4,1,2,0,6,4,13,3,1,4,0}, 
    {1,1,7,0,7,4,6,2,5,5,4,6,0,5,5,9,2,5,3,3,5,5,2,6,9,0,2,4,9,4,7,5}, 
    {9,6,2,4,1,2,2,2,5,2,1,2,7,6,1,2,8,1,6,3,1,7,0,0,3,9,3,9,0,4,0,3}, 
    {4,0,3,3,6,3,3,4,0,3,5,8,0,8,6,4,7,2,0,6,3,2,1,4,5,0,4,7,9,2,5,4}, 
    {10,6,0,4,5,5,2,2,3,0,2,3,5,9,0,5,5,0,4,6,3,3,0,3,4,8,2,6,3,6,0,7}, 
    {4,4,2,4,3,6,6,2,0,3,5,7,2,6,3,7,4,0,1,5,1,0,0,5,0,5,3,11,7,4,3,3}, 
    {1,5,5,1,5,7,9,0,4,5,3,5,3,3,3,12,0,2,4,3,3,3,2,7,4,5,0,8,8,6,5,4}, 
    {0,4,7,0,1,5,9,0,6,7,2,3,4,0,3,9,2,3,7,0,1,7,2,4,3,5,2,11,4,4,6,0}
};

#ifdef FOR_BACK
int16_t result[16][32] = {
    {9,5,9,4,7,3,4,5,13,6,4,3,8,6,4,6,7,9,11,5,7,13,6,3,13,5,4,5,5,5,6,7},
    {2,2,11,1,7,4,7,6,9,7,7,8,4,1,8,6,6,7,7,1,5,9,6,6,10,2,8,8,7,3,9,3},
    {10,3,6,6,4,1,1,6,8,4,3,7,7,10,4,1,12,5,7,8,6,11,3,3,10,8,8,5,4,3,4,5},
    {10,5,7,0,11,2,4,4,12,5,5,4,6,8,3,9,4,6,7,7,10,10,4,6,14,6,5,2,8,9,5,10},
    {8,8,7,3,6,4,7,3,12,6,2,0,10,4,1,7,4,6,11,3,5,11,1,5,8,9,5,8,0,7,4,5},
    {6,3,7,6,3,1,4,8,3,6,5,8,3,7,9,3,10,5,4,6,3,8,3,3,5,4,7,10,7,0,7,2},
    {12,2,3,6,8,3,1,7,6,1,5,5,5,11,5,2,9,5,3,9,6,6,3,3,11,3,4,5,7,4,3,10},
    {12,9,5,2,9,7,7,2,9,2,2,4,8,8,2,10,2,4,9,5,8,8,5,8,9,12,1,5,7,11,3,9},
    {5,6,6,5,1,6,6,6,4,6,8,6,3,4,5,4,9,2,4,4,2,5,5,3,1,9,7,14,6,1,6,1},
    {5,3,8,4,9,6,7,3,8,7,8,8,2,5,6,9,4,7,5,7,6,6,4,8,10,3,4,7,12,6,8,6},
    {12,9,5,7,3,3,3,5,8,5,5,5,9,6,3,4,8,1,8,7,5,8,4,1,5,12,4,12,3,4,2,4},
    {8,3,5,6,8,5,3,5,2,6,8,11,0,9,7,6,9,3,1,10,6,0,5,5,6,1,7,9,11,4,7,5},
    {14,9,0,6,5,7,2,4,5,1,6,5,7,9,1,7,7,2,6,8,6,3,0,5,6,10,5,8,5,7,1,9},
    {7,6,4,7,5,7,7,4,2,5,9,9,4,6,5,8,5,1,3,9,4,2,1,8,1,7,5,14,9,3,5,5},
    {1,8,6,4,7,6,10,3,6,7,6,8,5,4,4,14,1,4,6,5,5,4,6,9,6,9,4,9,10,8,7,5},
    {1,7,9,2,4,6,10,2,8,10,3,7,5,2,4,11,4,4,9,1,4,9,3,5,4,6,3,14,7,4,6,1}
};
#else
int16_t result[16][32] = {
    {12,7,11,6,9,5,4,5,17,6,4,4,9,4,6,4,6,10,12,6,8,13,5,2,14,6,3,5,3,4,8,7},
    {3,5,10,1,5,8,6,6,11,7,6,6,5,0,8,6,5,13,6,0,4,10,8,5,12,1,10,9,8,4,9,2},
    {14,2,5,6,3,2,1,5,10,2,3,6,5,10,4,0,17,8,6,7,5,11,4,5,8,7,7,6,3,3,4,4},
    {8,5,7,0,10,2,3,6,12,4,5,3,10,8,4,10,3,6,7,7,11,9,7,5,13,9,3,1,10,9,5,9},
    {9,9,6,4,4,6,6,2,11,7,4,1,8,4,0,7,3,6,11,4,4,13,1,5,7,9,3,11,1,10,4,6},
    {6,1,7,7,3,0,6,7,6,5,9,8,4,6,12,4,9,6,6,3,3,7,3,3,4,7,6,10,10,0,7,1},
    {11,5,3,6,8,4,2,10,5,3,4,4,3,12,6,2,8,4,2,10,6,7,2,2,9,3,3,6,6,4,3,8},
    {12,14,4,6,10,9,6,5,8,2,4,4,9,8,4,12,6,6,8,4,7,7,6,10,8,12,0,4,8,9,2,8},
    {2,9,7,8,0,5,6,4,5,8,8,6,3,2,4,4,9,2,4,2,1,8,3,3,3,11,6,14,5,5,5,0},
    {4,4,7,9,8,9,7,3,9,12,8,8,3,6,6,11,6,7,4,8,5,8,5,9,9,3,2,5,13,4,8,5},
    {9,7,3,9,8,4,6,7,10,7,9,3,7,6,6,5,8,1,7,9,7,8,6,0,3,15,3,10,3,4,5,3},
    {7,5,4,6,7,8,4,6,1,7,8,9,0,10,9,8,7,2,3,9,4,1,4,5,11,3,5,13,10,6,7,5},
    {16,10,1,7,6,9,2,3,6,0,7,6,5,8,3,7,5,4,5,6,4,3,0,4,4,12,4,9,4,10,0,7},
    {8,8,3,10,6,10,11,3,5,3,8,8,2,5,8,7,4,3,4,12,5,4,3,13,0,6,3,14,8,3,8,4},
    {1,7,10,6,6,6,9,7,7,8,4,6,3,2,3,15,0,6,5,5,4,4,6,8,8,11,2,9,9,11,7,4},
    {4,5,11,4,5,7,13,1,7,10,2,7,6,2,10,13,4,3,11,0,5,8,8,7,3,9,2,15,7,3,6,2}
};
#endif

#ifdef FOR_BACK
    #ifndef VECTOR_EXT
int16_t wires[q_field][q_field] = {
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
    {1,0,5,9,15,2,11,14,10,3,8,6,13,12,7,4},
    {2,5,0,6,10,1,3,12,15,11,4,9,7,14,13,8},
    {3,9,6,0,7,11,2,4,13,1,12,5,10,8,15,14},
    {4,15,10,7,0,8,12,3,5,14,2,13,6,11,9,1},
    {5,2,1,11,8,0,9,13,4,6,15,3,14,7,12,10},
    {6,11,3,2,12,9,0,10,14,5,7,1,4,15,8,13},
    {7,14,12,4,3,13,10,0,11,15,6,8,2,5,1,9},
    {8,10,15,13,5,4,14,11,0,12,1,7,9,3,6,2},
    {9,3,11,1,14,6,5,15,12,0,13,2,8,10,4,7},
    {10,8,4,12,2,15,7,6,1,13,0,14,3,9,11,5},
    {11,6,9,5,13,3,1,8,7,2,14,0,15,4,10,12},
    {12,13,7,10,6,14,4,2,9,8,3,15,0,1,5,11},
    {13,12,14,8,11,7,15,5,3,10,9,4,1,0,2,6},
    {14,7,13,15,9,12,8,1,6,4,11,10,5,2,0,3},
    {15,4,8,14,1,10,13,9,2,7,5,12,11,6,3,0}
};
    #else
// Index = wires * 2 bytes * 4 columns (dc)
int16_t wires_idx[q_field][q_field] = {
    {0,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120},
    {8,0,40,72,120,16,88,112,80,24,64,48,104,96,56,32},
    {16,40,0,48,80,8,24,96,120,88,32,72,56,112,104,64},
    {24,72,48,0,56,88,16,32,104,8,96,40,80,64,120,112},
    {32,120,80,56,0,64,96,24,40,112,16,104,48,88,72,8},
    {40,16,8,88,64,0,72,104,32,48,120,24,112,56,96,80},
    {48,88,24,16,96,72,0,80,112,40,56,8,32,120,64,104},
    {56,112,96,32,24,104,80,0,88,120,48,64,16,40,8,72},
    {64,80,120,104,40,32,112,88,0,96,8,56,72,24,48,16},
    {72,24,88,8,112,48,40,120,96,0,104,16,64,80,32,56},
    {80,64,32,96,16,120,56,48,8,104,0,112,24,72,88,40},
    {88,48,72,40,104,24,8,64,56,16,112,0,120,32,80,96},
    {96,104,56,80,48,112,32,16,72,64,24,120,0,8,40,88},
    {104,96,112,64,88,56,120,40,24,80,72,32,8,0,16,48},
    {112,56,104,120,72,96,64,8,48,32,88,80,40,16,0,24},
    {120,32,64,112,8,80,104,72,16,56,40,96,88,48,24,0}
};
    #endif
#endif

int main(void)
{   
    int i,j,k,var,row; /* Used for loop indexs */

    int16_t 	decoded[N_code][m_field];
    int16_t 	Rmn_SRL[M_code][q_field][dc];
    int16_t Qmn[q_field][dc];
    int16_t aux1, aux2;
    int16_t MAX_temp;          
    int16_t Rmn[q_field][dc];
    


    int16_t min_vec[q_field];

    int16_t Qn_NEW[q_field][dc];   
    int16_t Qmn_temp[q_field]; 


#ifndef FOR_BACK
    int16_t min_temp;
	int16_t min_temp2;
    int16_t min1[q_field], min2[q_field];
    int16_t cam1[q_field], cam2[q_field];
    int16_t dQ_min1, dQ_min2; 
    int16_t dQ_pos1, dQ_pos2;
    int16_t dWmn[q_field][dc], dWmn2[q_field][dc];
    int16_t temp[q_field][dc], temp2[q_field];
    int16_t pos[q_field];
    int16_t cam1_temp[q_field/2 - 1], cam2_temp[q_field/2 - 1];
    int16_t min_global[q_field]; //min2_global[q_field];  
    int16_t oRmn_temp;
    int16_t pos_temp;
    int16_t max[q_field/2];
    int16_t beta;
    int16_t z[dc];
#endif    
	
#ifdef FOR_BACK
	int a;
	int16_t R_Forward[q_field][dc], R_Backward[q_field][dc];
    
    

    #ifdef VECTOR_EXT
    int16_t R_compare_vec[q_field];
    #else
    int16_t R_compare, R_compare_2, R_compare_3;
    int16_t min_temp, min_temp_2, min_temp_3;
    int wires_aux;
    #endif
#endif    

    printf("--------------------------\n");
    printf("    NB LDPC Algorithm\n");
    printf("--------------------------\n");
    unsigned long int block_size_p;
    // Set the vector configuration
    asm volatile("vsetvli %0, %1, e16, m1, ta, ma" : "=r"(block_size_p) : "r"(q_field));

    /* Inicialize Rmn_SRL to zero values */
    for (i = 0; i<M_code; i++)
        for (j=0; j<q_field; j++)
            for (k=0; k<dc; k++)
                Rmn_SRL[i][j][k] = 0;
	
    
    /* Loop through all decoding iterations */
    for (var=0; var<it_max; var++)
    {						
        /* Loop through all rows of parity check matrix */
        for (row=0; row<M_code; row++)
        {
#ifndef FOR_BACK            
            beta = 0; /* Inicializa el valor del sindrome a cero */
#endif            
            /* Extract Qmn(a) messages from Qn(a) memories */
            for (i=0;i<dc;i++)
            {
                aux2 = pow_coefH[row][i] + 1;
                aux1 = q_field - aux2;

                /* PERMUTACION DE LA COLUMNA EXTRAIDA */

#ifndef VECTOR_EXT                
                for (j=0;j<q_field;j++)
                    Qmn_temp[j] = Ln_aux[j][col[row][i]];            


                if (pow_coefH[row][i]==0)
                {
                    for (j=0;j<q_field;j++)
                        Qmn[j][i] = Qmn_temp[j];                 
                }
                else
                {		                
                    for (k = aux2; k < q_field; k++)
                    {
                        Qmn[k][i] = Qmn_temp[k-aux2+1];
                    }
                    for (k = 1; k < aux2; k++)
                    {
                        Qmn[k][i] = Qmn_temp[k+aux1];
                    }
                    Qmn[0][i] = Qmn_temp[0];
                }                              

                /***********************************/

                /* Qmn_prima = Qmn_p - reshape(Rmn_SRL(i,:,:),q_field,dc); */

                MAX_temp = 10000;
                for (j=0;j<q_field;j++)
                {
                    Qmn[j][i] = Qmn[j][i] - Rmn_SRL[row][j][i];

                    /* Busqueda del minimo para normalizacion */
                    if (Qmn[j][i] < MAX_temp) {
                        MAX_temp = Qmn[j][i];
    #ifndef FOR_BACK                        
                        z[i] = pos_exp_table[j];  
    #endif                        
                    }
                }

    #ifndef FOR_BACK                   
                beta = gfadd[beta][z[i]];  /* SINDROME DEL CHECK NODE */
    #endif
                /*****************************************/

                /* Normalizacion y precision finita*/
                for (j=0;j<q_field;j++)
                    Qmn[j][i] = Qmn[j][i]-MAX_temp; 

                /*****************/
#endif
#ifdef VECTOR_EXT
                /***************************************
                ****************************************
                **************** PERMUTE ***************
                ****************************************
                ***************************************/
                // 16 (int16_t) * 32 (columns) / 8 
                asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&Ln_aux[0][col[row][i]]), "r"(64));
                if (pow_coefH[row][i]!=0)
                {
                    // Clean vectors     
                    asm volatile("vmv.v.i v1, 0;");
                    asm volatile("vmv.v.i v2, 0;");
                    asm volatile("vmv.v.i v3, 0;");
                    asm volatile("vmv.v.i v4, 0;");
                    asm volatile("vmv.v.i v5, 0;");
                    asm volatile("vmv.v.i v6, 0;");
                    asm volatile("vmv.v.i v7, 0;");

                    // Permutation and clean position 0
                    asm volatile("vslidedown.vi v1, v0, 1;");
                    asm volatile("vslideup.vx v2, v1, %0;" :: "r"(aux2));                
                    asm volatile("vslidedown.vx v3, v0, %0;" :: "r"(aux1+1));
                    asm volatile("vslideup.vx v4, v3, %0;" :: "r"(1));   
                    asm volatile("vadd.vv v5, v2, v4;");

                    // Set position 0
                    asm volatile("vmv.v.x v6, %0;" :: "r"(Ln_aux[0][col[row][i]]));
                    asm volatile("vslidedown.vx v7, v6, %0;" :: "r"(q_field-1));

                    // Merge both vectors
                    asm volatile("vadd.vv v0, v5, v7;");
                }                 
                asm volatile("vsse16.v v0, (%0), %1;" ::"r"(&Qmn[0][i]), "r"(8)); 

                /***************************************
                ****************************************
                **************** MINIMUM ***************
                ****************************************
                ***************************************/
                // OFFSET
                // 16 (int16_t) * 4 (columns) / 8 
                asm volatile("vlse16.v v12, (%0), %1;" ::"r"(&Rmn_SRL[row][0][i]), "r"(8));                 
                asm volatile("vsub.vv v8, v0, v12;");
                      
        #ifdef FOR_BACK  
                // Find minimum
                asm volatile("vmv.v.i v9, 15;");
                asm volatile("vredmin.vs v10, v8, v9;"); 
                asm volatile("vse16.v v10, (%0);" ::"r"(&min_vec[0])); 
                // SUBSTRACT OFFSET
                asm volatile("vsub.vx v11, v8, %0;" :: "r"(min_vec[0]));

        #else 
                //asm volatile("vsse16.v v8, (%0), %1;" ::"r"(&Qmn[0][i]), "r"(8));   
                MAX_temp = 10000;
                for (j=0;j<q_field;j++)
                {
                    /* Busqueda del minimo para normalizacion */
                    if (Qmn[j][i] < MAX_temp) {
                        MAX_temp = Qmn[j][i];                      
                        z[i] = pos_exp_table[j];                      
                    }
                }
                beta = gfadd[beta][z[i]];  /* SINDROME DEL CHECK NODE */
                // SUBSTRACT OFFSET
                asm volatile("vsub.vx v11, v8, %0;" :: "r"(MAX_temp));
        #endif
                // 16 (int16_t) * 4 (columns) / 8 
                asm volatile("vsse16.v v11, (%0), %1;" ::"r"(&Qmn[0][i]), "r"(8));                          
#endif  
            } 

            /*for(j=0;j<dc;j++){
                for(i=0;i<q_field;i++){
                    if(Qmn_aux[i][j] !=Qmn[i][j])
                        printf("ERROR\n");
                }
            }*/

            /*******************************************/
            /**********  CNU FUNCTION ******************/
            /*******************************************/
#ifdef FOR_BACK
            // Clean Forward and Backward variables and set value to first and last column
            // Clean R_Aux   
        #ifndef VECTOR_EXT    
            for(i=0;i<q_field;i++){        
                R_Forward[i][0] = Qmn[i][0];
                R_Backward[i][dc-1] = Qmn[i][dc-1];
                Rmn[i][dc-1] = Qmn[i][dc-1];
            }    
        #else
            // 16 (int16_t) * 4 (columns) / 8

            asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&Qmn[0][0]), "r"(8)); 
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Forward[0][0]), "r"(8));

            asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&Qmn[0][dc-1]), "r"(8)); 
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Backward[0][dc-1]), "r"(8));
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Rmn[0][dc-1]), "r"(8));
        #endif 

            for(a=2;a<dc;a++){
        #ifdef VECTOR_EXT
                asm volatile("vmv.v.i v25, 15;"); 

                // 16 (int16_t) * 16 (columns) / 8
                asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&wires_idx[0][0]), "r"(32));          
                asm volatile("vlse16.v v1, (%0), %1;" ::"r"(&wires_idx[0][1]), "r"(32));
                asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&wires_idx[0][2]), "r"(32));
                asm volatile("vlse16.v v3, (%0), %1;" ::"r"(&wires_idx[0][3]), "r"(32));
                asm volatile("vlse16.v v4, (%0), %1;" ::"r"(&wires_idx[0][4]), "r"(32));
                asm volatile("vlse16.v v5, (%0), %1;" ::"r"(&wires_idx[0][5]), "r"(32));
                asm volatile("vlse16.v v6, (%0), %1;" ::"r"(&wires_idx[0][6]), "r"(32));
                asm volatile("vlse16.v v7, (%0), %1;" ::"r"(&wires_idx[0][7]), "r"(32));  

                /***************************************
                ****************************************
                ***************** FORWARD **************
                ****************************************
                ***************************************/
                
                // Load Qmn[0][a-1]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&Qmn[0][a-1])); 
                
                // Load R_Forward[k][a-2]
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Forward[0][a-2]), "r"(8));
                
                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");
                
                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 
                
                // Store 
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[0][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[1][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[2][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[3][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[4][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[5][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[6][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[7][a-1] = R_compare_vec[0];

                /***************************************
                ****************************************
                **************** BACKWARD **************
                ****************************************
                ***************************************/
                
                // Load Qmn[wires_aux][4-a]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&Qmn[0][4-a])); 
                
                // Load R_Backward[k][5-a]
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Backward[0][5-a]), "r"(8));
                
                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");
                
                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 
                
                // Store 
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[0][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[1][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[2][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[3][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[4][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[5][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[6][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[7][4-a] = R_compare_vec[0];

                /***************************************
                ****************************************
                *************** PART 2 *****************
                ****************************************
                ***************************************/

                // 16 (int16_t) * 16 (columns) / 8
                asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&wires_idx[0][8]), "r"(32));          
                asm volatile("vlse16.v v1, (%0), %1;" ::"r"(&wires_idx[0][9]), "r"(32));
                asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&wires_idx[0][10]), "r"(32));
                asm volatile("vlse16.v v3, (%0), %1;" ::"r"(&wires_idx[0][11]), "r"(32));
                asm volatile("vlse16.v v4, (%0), %1;" ::"r"(&wires_idx[0][12]), "r"(32));
                asm volatile("vlse16.v v5, (%0), %1;" ::"r"(&wires_idx[0][13]), "r"(32));
                asm volatile("vlse16.v v6, (%0), %1;" ::"r"(&wires_idx[0][14]), "r"(32));
                asm volatile("vlse16.v v7, (%0), %1;" ::"r"(&wires_idx[0][15]), "r"(32)); 

                /***************************************
                ****************************************
                ***************** FORWARD **************
                ****************************************
                ***************************************/

                // Load Qmn[0][a-1]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&Qmn[0][a-1])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&Qmn[0][a-1])); 
                
                // Load R_Forward[k][a-2]
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Forward[0][a-2]), "r"(8));

                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");

                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 
                
                // Store 
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[8][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[9][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[10][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[11][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[12][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[13][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[14][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Forward[15][a-1] = R_compare_vec[0];

                /***************************************
                ****************************************
                **************** BACKWARD **************
                ****************************************
                ***************************************/
                
                // Load Qmn[wires_aux][4-a]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&Qmn[0][4-a])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&Qmn[0][4-a])); 
                
                // Load R_Backward[k][5-a]
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Backward[0][5-a]), "r"(8));

                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");

                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 
                
                // Store 
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[8][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[9][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[10][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[11][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[12][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[13][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[14][4-a] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                R_Backward[15][4-a] = R_compare_vec[0];

        #else
                for(i=0;i<q_field;i++){
                    min_temp = 10000;
                    min_temp_2 = 10000;

                    for(k=0;k<q_field;k++){
                        // Search maximum
                        wires_aux = wires[k][i];                
                        // Forward
                        if(R_Forward[k][a-2] > Qmn[wires_aux][a-1])
                            R_compare = R_Forward[k][a-2];                            
                        else
                            R_compare = Qmn[wires_aux][a-1];

                        // Backward
                        if(R_Backward[k][5-a] > Qmn[wires_aux][4-a])
                            R_compare_2 = R_Backward[k][5-a];
                        else
                            R_compare_2 = Qmn[wires_aux][4-a];

                        // Search minimum
                        // Forward
                        if(R_compare < min_temp)
                            min_temp = R_compare;
                        // Backward
                        if(R_compare_2 < min_temp_2)
                            min_temp_2 = R_compare_2;
                    }
                    R_Forward[i][a-1] = min_temp;
                    R_Backward[i][4-a] = min_temp_2;
                    //printf("R_Forward[%d][%d] = %d\n",i, a-1, min_temp);   
                }
        #endif

        #ifdef VECTOR_EXT
                
                // Load wire index in v0
                // 16 (int16_t) * 16 (columns) / 8
                /*asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&wires_idx[0][8]), "r"(32));          
                asm volatile("vlse16.v v1, (%0), %1;" ::"r"(&wires_idx[0][9]), "r"(32));
                asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&wires_idx[0][10]), "r"(32));
                asm volatile("vlse16.v v3, (%0), %1;" ::"r"(&wires_idx[0][11]), "r"(32));
                asm volatile("vlse16.v v4, (%0), %1;" ::"r"(&wires_idx[0][12]), "r"(32));
                asm volatile("vlse16.v v5, (%0), %1;" ::"r"(&wires_idx[0][13]), "r"(32));
                asm volatile("vlse16.v v6, (%0), %1;" ::"r"(&wires_idx[0][14]), "r"(32));
                asm volatile("vlse16.v v7, (%0), %1;" ::"r"(&wires_idx[0][15]), "r"(32)); */

                // Load R_Backward[wires_aux][a]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R_Backward[0][a])); 

                // Load R_Forward[k][a-2]
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Forward[0][a-2]), "r"(8));

                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");

                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 

                // Store 
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[8][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[9][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[10][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[11][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[12][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[13][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[14][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[15][a-1] = R_compare_vec[0];


                // 16 (int16_t) * 16 (columns) / 8
                asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&wires_idx[0][0]), "r"(32));          
                asm volatile("vlse16.v v1, (%0), %1;" ::"r"(&wires_idx[0][1]), "r"(32));
                asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&wires_idx[0][2]), "r"(32));
                asm volatile("vlse16.v v3, (%0), %1;" ::"r"(&wires_idx[0][3]), "r"(32));
                asm volatile("vlse16.v v4, (%0), %1;" ::"r"(&wires_idx[0][4]), "r"(32));
                asm volatile("vlse16.v v5, (%0), %1;" ::"r"(&wires_idx[0][5]), "r"(32));
                asm volatile("vlse16.v v6, (%0), %1;" ::"r"(&wires_idx[0][6]), "r"(32));
                asm volatile("vlse16.v v7, (%0), %1;" ::"r"(&wires_idx[0][7]), "r"(32));   

                // Load R_Backward[wires_aux][a]
                asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R_Backward[0][a])); 
                asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R_Backward[0][a])); 

                // Load R_Forward[k][a-2]
                // 16 (int16_t) * 4 (columns) / 8
                //asm volatile("vlse16.v v16, (%0), %1;" ::"r"(&R_Forward[0][a-2]), "r"(8));

                // Max find
                asm volatile("vmax.vv v17, v16, v8;");
                asm volatile("vmax.vv v18, v16, v9;");
                asm volatile("vmax.vv v19, v16, v10;");
                asm volatile("vmax.vv v20, v16, v11;");
                asm volatile("vmax.vv v21, v16, v12;");
                asm volatile("vmax.vv v22, v16, v13;");
                asm volatile("vmax.vv v23, v16, v14;");
                asm volatile("vmax.vv v24, v16, v15;");

                // Min find
                asm volatile("vredmin.vs v17, v17, v25;"); 
                asm volatile("vredmin.vs v18, v18, v25;"); 
                asm volatile("vredmin.vs v19, v19, v25;"); 
                asm volatile("vredmin.vs v20, v20, v25;"); 
                asm volatile("vredmin.vs v21, v21, v25;"); 
                asm volatile("vredmin.vs v22, v22, v25;"); 
                asm volatile("vredmin.vs v23, v23, v25;"); 
                asm volatile("vredmin.vs v24, v24, v25;"); 

                // Store 
                // 16 (int16_t) * 4 (columns) / 8
                asm volatile("vsse16.v v17, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[0][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v18, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[1][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v19, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[2][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v20, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[3][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[4][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[5][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[6][a-1] = R_compare_vec[0];
                asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&R_compare_vec[0]), "r"(8));
                Rmn[7][a-1] = R_compare_vec[0];

        #else
                // Standard Min Max
                for(i=0;i<q_field;i++){
                    min_temp_3 = 10000;        
                    // Search maximum
                    for(k=0;k<q_field;k++){
                        wires_aux = wires[k][i];
                        if(R_Forward[k][a-2] > R_Backward[wires_aux][a])
                                R_compare_3 = R_Forward[k][a-2];
                            else
                                R_compare_3 = R_Backward[wires_aux][a];
                        // Search minimum
                        if(R_compare_3 < min_temp_3)
                            min_temp_3 = R_compare_3;                
                    }        
                    Rmn[i][a-1] = min_temp_3;     
                    //printf("Rmn[%d][%d] = %d\n",i, a-1, min_temp_3);   
                }      
        #endif        
                    
            }

        #ifndef VECTOR_EXT
            for(i=0;i<q_field;i++){
                Rmn[i][0] =       R_Backward[i][1];
                Rmn[i][dc-1] =    R_Forward[i][dc-2];
            }
        #else
            // 16 (int16_t) * 4 (columns) / 8
            asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R_Backward[0][1]), "r"(8)); 
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Rmn[0][0]), "r"(8));

            asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R_Forward[0][dc-2]), "r"(8)); 
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Rmn[0][dc-1]), "r"(8));
        #endif

        /*if(row == 1){
            for(j=0;j<dc;j++){
                for(i=0;i<q_field;i++){
                    printf("%d, ", Rmn[i][j]);
                }
                printf("\n");
            }
        }*/

#else                                                        
            /* paso de mensajes al dominio delta */
        #ifdef VECTOR_EXT

            asm volatile("vle16.v v0, (%0);" ::"r"(pos_exp_table_idx));  
            asm volatile("vluxei16.v v1, (%0), v0;" ::"r"(&gfadd[z[0]][0])); 
            asm volatile("vluxei16.v v2, (%0), v0;" ::"r"(&gfadd[z[1]][0])); 
            asm volatile("vluxei16.v v3, (%0), v0;" ::"r"(&gfadd[z[2]][0])); 
            asm volatile("vluxei16.v v4, (%0), v0;" ::"r"(&gfadd[z[3]][0])); 
            
            // 16 (int16_t) * 4 (columns) / 8
            asm volatile("vsse16.v v1, (%0), %1;" ::"r"(&temp[0][0]), "r"(8)); 
            asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&temp[0][1]), "r"(8)); 
            asm volatile("vsse16.v v3, (%0), %1;" ::"r"(&temp[0][2]), "r"(8)); 
            asm volatile("vsse16.v v4, (%0), %1;" ::"r"(&temp[0][3]), "r"(8)); 
            
            asm volatile("vmul.vx v5, v1, %0" :: "r"(2)); 
            asm volatile("vmul.vx v6, v2, %0" :: "r"(2)); 
            asm volatile("vmul.vx v7, v3, %0" :: "r"(2)); 
            asm volatile("vmul.vx v8, v4, %0" :: "r"(2)); 

            asm volatile("vluxei16.v v9, (%0), v5;" ::"r"(exp_pos_table)); 
            asm volatile("vluxei16.v v10, (%0), v6;" ::"r"(exp_pos_table)); 
            asm volatile("vluxei16.v v11, (%0), v7;" ::"r"(exp_pos_table)); 
            asm volatile("vluxei16.v v12, (%0), v8;" ::"r"(exp_pos_table)); 

            asm volatile("vsub.vx v13, v9, %0" :: "r"(1)); 
            asm volatile("vsub.vx v14, v10, %0" :: "r"(1)); 
            asm volatile("vsub.vx v15, v11, %0" :: "r"(1)); 
            asm volatile("vsub.vx v16, v12, %0" :: "r"(1)); 

            asm volatile("vmul.vx v17, v13, %0" :: "r"(8)); 
            asm volatile("vmul.vx v18, v14, %0" :: "r"(8)); 
            asm volatile("vmul.vx v19, v15, %0" :: "r"(8)); 
            asm volatile("vmul.vx v20, v16, %0" :: "r"(8)); 

            asm volatile("vluxei16.v v21, (%0), v17;" ::"r"(&Qmn[0][0])); 
            asm volatile("vluxei16.v v22, (%0), v18;" ::"r"(&Qmn[0][1])); 
            asm volatile("vluxei16.v v23, (%0), v19;" ::"r"(&Qmn[0][2])); 
            asm volatile("vluxei16.v v24, (%0), v20;" ::"r"(&Qmn[0][3])); 

            // 16 (int16_t) * 4 (columns) / 8
            asm volatile("vsse16.v v21, (%0), %1;" ::"r"(&dWmn[0][0]), "r"(8)); 
            asm volatile("vsse16.v v22, (%0), %1;" ::"r"(&dWmn[0][1]), "r"(8)); 
            asm volatile("vsse16.v v23, (%0), %1;" ::"r"(&dWmn[0][2]), "r"(8)); 
            asm volatile("vsse16.v v24, (%0), %1;" ::"r"(&dWmn[0][3]), "r"(8)); 

        #else
            for (i=0;i<dc;i++)
            {
                for (j=0;j<q_field;j++)
                {
                    temp[j][i] = gfadd[z[i]][pos_exp_table[j]];
                    dWmn[j][i] = Qmn[ exp_pos_table[ temp[j][i] ] -1][i];
                }             

            }
        #endif

            /* Busqueda de minimos */
            for (j=0;j<q_field;j++)
            {
                min_temp = dWmn[j][0];
                min_temp2 = 10000;
                pos_temp = 0;
            
                for (i=1;i<dc;i++)
                {
                    if (dWmn[j][i] < min_temp)
                    {
                        min_temp2 = min_temp;
                        min_temp = dWmn[j][i];
                        pos_temp = i;
                    }
                    else if (dWmn[j][i] < min_temp2)
                    {
                        min_temp2 = dWmn[j][i];
                    }
                }

                min1[j] = min_temp;
                min2[j] = min_temp2;
                pos[j] = pos_temp;
            }   

            /* Calculo de la columna extra */
            for (i=0;i<q_field-1;i++)
            {
                min_temp = min1[i+1];
                //float min2_temp = 10000;
                pos_temp = 0;
                
                max[0] = min1[i+1];
                for (j=0;j<q_field/2-1;j++)
                {
                    cam1_temp[j] = pos[conf_tb[j][0][i]-1];
                    cam2_temp[j] = pos[conf_tb[j][1][i]-1];
                    if (cam1_temp[j] == cam2_temp[j]) 
                    {
                        max[j+1] = 10000;
                    }
                    else
                    {
                        //max[j+1] = min1[conf_tb[j][0][i]-1] + min1[conf_tb[j][1][i]-1];
                        
                        if (min1[conf_tb[j][0][i]-1] > min1[conf_tb[j][1][i]-1])
                        {
                            max[j+1] = min1[conf_tb[j][0][i]-1];
                        }
                        else
                        {
                            max[j+1] = min1[conf_tb[j][1][i]-1];
                        }
                        
                    }
                    
                    if (max[j+1]<min_temp)
                    {
                        min_temp = max[j+1];
                        pos_temp = j+1;
                    }
                }
                
                min_global[i+1] = min_temp;
                
                if (pos_temp == 0)
                {
                    cam1[i+1] = pos[i+1];
                    cam2[i+1] = -1;
                }
                else
                {
                    cam1[i+1] = cam1_temp[pos_temp-1];
                    cam2[i+1] = cam2_temp[pos_temp-1];
                }
            }
                                        
            min_global[0] = 0;
                  
            /* Busqueda de los dos minimos valores en la columna extra */
            dQ_min1 = 10000;
            dQ_min2 = 10000;
            dQ_pos1 = 0;
            dQ_pos2 = 0;
            for (j=1;j<q_field;j++)
            {
                if (min_global[j] < dQ_min1)
                {
                    dQ_min2 = dQ_min1;
                    dQ_min1 = min_global[j];
                    dQ_pos2 = dQ_pos1;
                    dQ_pos1 = j;
                }
                else if (min_global[j] < dQ_min2)
                {
                    dQ_min2 = min_global[j];
                    dQ_pos2 = j;
                }
            }

#endif            
            for (i=0;i<dc;i++)
            {				      
#ifndef FOR_BACK                          
                for (j=0;j<q_field;j++)
                {
                    if ( (j == dQ_pos1) || (j == dQ_pos2) || (j == 0) ) /* Si estamos en el elemento del campo que contiene el minimo de la columna extra entonces... */
                    {
                        if ((cam2[j] == -1)&&(i == cam1[j]))
                            oRmn_temp = min2[j];
                        else if ((i == cam1[j]) || (i == cam2[j]))
                            oRmn_temp = min1[j];
                        else
                            oRmn_temp = min_global[j];

                    }
                    else
                    {
                        //oRmn_temp = voto * MAX_temp2;
                        if ((cam2[j] == -1)&&(i == cam1[j]))
                            oRmn_temp = min2[j];
                        else if ((i == cam1[j]) || (i == cam2[j]))
                            oRmn_temp = min1[j];
                        else
                            oRmn_temp = dQ_min2;// voto * dQ_min2;
                    }
                    
                    temp[j][i] = exp_pos_table[gfadd[temp[j][i]][beta]]-1;
                    
                    /*Aplica escalado y precision finita a los mensajes de salida del CN*/
                    Rmn[temp[j][i]][i] = oRmn_temp;//*scaling_factor;
                }
#endif

#ifdef VECTOR_EXT  

                // 16 (int16_t) * 4 (columns) / 8 
                asm volatile("vlse16.v v0, (%0), %1;" ::"r"(&Rmn[0][i]), "r"(8));
                asm volatile("vlse16.v v1, (%0), %1;" ::"r"(&Qmn[0][i]), "r"(8));
                asm volatile("vadd.vv v2, v0, v1;");
                asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Rmn[0][i]), "r"(8));
                asm volatile("vsse16.v v0, (%0), %1;" ::"r"(&Rmn_SRL[row][0][i]), "r"(8));
#else
                for (j=0;j<q_field;j++)
                {
                
                    Rmn_SRL[row][j][i] = Rmn[j][i];

                    /*CNout = CNout + Qmn_prima;*/
                    Rmn[j][i] = Rmn[j][i] + Qmn[j][i];
                }
#endif

                /****** PERMUTACION INVERSA *******/
                aux1 = pow_coefH[row][i];
                aux2 = q_field - aux1;
#ifdef VECTOR_EXT             
                asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&Rmn[0][i]), "r"(8));
                if (aux1!=0)
                {
                    // Clean vectors     
                    asm volatile("vmv.v.i v0, 0;");
                    asm volatile("vmv.v.i v4, 0;");
                    asm volatile("vmv.v.i v6, 0;");
                    asm volatile("vmv.v.i v8, 0;");
                    asm volatile("vmv.v.i v10, 0;");

                    // Permutation and clean position 0
                    asm volatile("vslidedown.vi v0, v2, 1;");
                    asm volatile("vslideup.vx v4, v0, %0;" :: "r"(aux2));                
                    asm volatile("vslidedown.vx v6, v2, %0;" :: "r"(aux1+1));
                    asm volatile("vslideup.vx v8, v6, %0;" :: "r"(1));   
                    asm volatile("vadd.vv v6, v4, v8;");

                    // Set position 0
                    asm volatile("vmv.v.i v8, 0;");
                    asm volatile("vmv.v.x v8, %0;" :: "r"(Rmn[0][i]));
                    asm volatile("vslidedown.vx v10, v8, %0;" :: "r"(q_field-1));

                    // Merge both vectors
                    asm volatile("vadd.vv v2, v6, v10;");
                }
                // 16 (int16_t) * 4 (columns) / 8 
                //asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Qn_NEW[0][i]), "r"(8));   

#else 

                if (aux1==0) /* Si hmn es diferente de 0 o 1 entonces se permuta*/
                {
                    for (j=0;j<q_field;j++)
                    {
                        Qn_NEW[j][i] = Rmn[j][i];
                    }
                }
                else
                {
                    for (k = 1; k < q_field-aux1; k++)
                    {
                        Qn_NEW[k][i] = Rmn[aux1+k][i];
                    }
                    for (k = q_field-aux1; k < q_field; k++)
                    {
                        Qn_NEW[k][i] = Rmn[k-q_field+aux1+1][i];
                    }

                    Qn_NEW[0][i] = Rmn[0][i];
                }
#endif                
 
                /***********************************/

                /* PRECISION FINITA Y ACTUALIZACION DE LAS VN*/
#ifdef VECTOR_EXT
                // 16 (int16_t) * 32 (columns) / 8 
                asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&Ln_aux[0][col[row][i]]), "r"(64));
#else
                for (j=0;j<q_field;j++)
                    Ln_aux[j][col[row][i]] = Qn_NEW[j][i];
#endif                                    
            }

        }

        int error = 0;
        for (i=0;i<16;i++){
            for (j=0;j<32;j++){
                //printf("%d,",Ln_aux[i][j]);
                if(Ln_aux[i][j]!= result[i][j])
                    error = 1;
            }
            //printf("\n");
        }


        if(error ==1)
            printf("ERROR\n");
        else 
            printf("SUCCESS\n");

        

        /*Tentatively decoding*/
        /*Compute: Qn(a)=Ln(a)+sum(Rmn)*/

        aux2 = 0;
        for (i=0;i<N_code;i++)
        {
            MAX_temp = 10000;
            for (j=0;j<q_field;j++)
            {
                /* Busqueda del minimo para normalizacion */
                if (Ln_aux[j][i] < MAX_temp)
                {
                    MAX_temp = Ln_aux[j][i];
                    aux1 = j;
                }
            }

            for (j=0;j<m_field;j++)
            {
                decoded[i][j] = fld_bpsk_bip[pos_exp_table[aux1]][j];
                if (decoded[i][j] != m[i][j])
                {
                    aux2 = aux2+1;
                }
            }
        }

        if (aux2 == 0)
            break;
    }

    asm volatile("jr x0;");
    return 0;
}