#include <stdint.h>
#include <runtime.h>
#include <printf.h>

#define q_field 16
#define m_field 4
#define N_code 32
#define M_code 16
#define dc 4
#define dv 2
#define VECTOR_EXT

int16_t R[q_field][dc] = {
    {6,11,6,12},
    {2,2,5,5},
    {1,4,10,4},
    {9,8,3,0},
    {4,1,3,4},
    {10,5,9,3},
    {4,5,8,8},
    {1,0,1,8},
    {0,3,7,13},
    {0,0,2,7},
    {9,4,8,4},
    {8,6,7,9},
    {6,7,5,8},
    {3,6,4,0},
    {11,9,0,9},
    {8,10,2,3}
};

int16_t result[q_field][dc] = {
    {1,1,1,1},
    {1,1,0,1},
    {2,0,1,1},
    {1,1,0,0},
    {2,1,0,1},
    {2,1,0,0},
    {1,1,0,1},
    {0,0,1,1},
    {1,1,1,1},
    {0,1,1,1},
    {1,1,0,0},
    {0,0,1,1},
    {0,1,1,1},
    {1,1,0,1},
    {1,1,0,0},
    {2,0,1,1}
};
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

int main(void) 
{
    printf("--------------------------\n");
    printf("Forward Backward algorithm\n");
    printf("--------------------------\n");
    
    int i, j, a;
    int16_t R_Forward[q_field][dc], R_Backward[q_field][dc], R_Aux[q_field][dc];

#ifndef VECTOR_EXT
    int16_t R_compare, R_compare_2, R_compare_3;
    int16_t min_temp, min_temp_2, min_temp_3, min3;
    int wires_aux;
    int k;
#else    
    int16_t R_compare_vec[q_field];
    unsigned long int block_size_p;
    // Vector configuration
    asm volatile("vsetvli %0, %1, e16, m1, ta, ma" : "=r"(block_size_p) : "r"(q_field));
#endif
    

    // Clean Forward and Backward variables and set value to first and last column
    // Clean R_Aux   
#ifndef VECTOR_EXT    
    for(i=0;i<q_field;i++){        
        R_Forward[i][0] = R[i][0];
        R_Backward[i][dc-1] = R[i][dc-1];
        R_Aux[i][dc-1] = R[i][dc-1];
    }    
#else
    // 16 (int16_t) * 4 (columns) / 8

    asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R[0][0]), "r"(8)); 
    asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Forward[0][0]), "r"(8));

    asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R[0][dc-1]), "r"(8)); 
    asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Backward[0][dc-1]), "r"(8));
    asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Aux[0][dc-1]), "r"(8));
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
         
        // Load R[0][a-1]
        asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R[0][a-1])); 
        
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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[0][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[1][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[2][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[3][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[4][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[5][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[6][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[7][a-1] = R_compare_vec[0];

        /***************************************
        ****************************************
        **************** BACKWARD **************
        ****************************************
        ***************************************/
         
        // Load R[wires_aux][4-a]
        asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R[0][4-a])); 
        
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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[0][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[1][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[2][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[3][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[4][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[5][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[6][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
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

        // Load R[0][a-1]
        asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R[0][a-1])); 
        asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R[0][a-1])); 
        
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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[8][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[9][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[10][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[11][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[12][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[13][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[14][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
        R_Forward[15][a-1] = R_compare_vec[0];

        /***************************************
        ****************************************
        **************** BACKWARD **************
        ****************************************
        ***************************************/
        
        // Load R[wires_aux][4-a]
        asm volatile("vluxei16.v v8, (%0), v0;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v9, (%0), v1;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v10, (%0), v2;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v11, (%0), v3;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v12, (%0), v4;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v13, (%0), v5;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v14, (%0), v6;" ::"r"(&R[0][4-a])); 
        asm volatile("vluxei16.v v15, (%0), v7;" ::"r"(&R[0][4-a])); 
        
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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[8][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[9][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[10][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[11][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[12][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[13][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[14][4-a] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
        R_Backward[15][4-a] = R_compare_vec[0];

#else
        for(i=0;i<q_field;i++){
            min_temp = 10000;
            min_temp_2 = 10000;

            for(k=0;k<q_field;k++){
                // Search maximum
                wires_aux = wires[k][i];                
                // Forward
                if(R_Forward[k][a-2] > R[wires_aux][a-1])
                    R_compare = R_Forward[k][a-2];                            
                else
                    R_compare = R[wires_aux][a-1];

                // Backward
                if(R_Backward[k][5-a] > R[wires_aux][4-a])
                    R_compare_2 = R_Backward[k][5-a];
                else
                    R_compare_2 = R[wires_aux][4-a];

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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[8][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[9][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[10][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[11][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[12][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[13][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[14][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[15][a-1] = R_compare_vec[0];


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
        asm volatile("vse16.v v17, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[0][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v18, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[1][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v19, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[2][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v20, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[3][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v21, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[4][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v22, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[5][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v23, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[6][a-1] = R_compare_vec[0];
        asm volatile("vse16.v v24, (%0);" ::"r"(&R_compare_vec[0]));
        R_Aux[7][a-1] = R_compare_vec[0];

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
            R_Aux[i][a-1] = min_temp_3;     
            //printf("R_Aux[%d][%d] = %d\n",i, a-1, min_temp_3);   
        }      
#endif        
               
    }

#ifndef VECTOR_EXT
    for(i=0;i<q_field;i++){
        R_Aux[i][0] =       R_Backward[i][1];
        R_Aux[i][dc-1] =    R_Forward[i][dc-2];
    }
#else
    // 16 (int16_t) * 4 (columns) / 8
    asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R_Backward[0][1]), "r"(8)); 
    asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Aux[0][0]), "r"(8));

    asm volatile("vlse16.v v2, (%0), %1;" ::"r"(&R_Forward[0][dc-2]), "r"(8)); 
    asm volatile("vsse16.v v2, (%0), %1;" ::"r"(&R_Aux[0][dc-1]), "r"(8));
#endif

    // Print
    int error = 0;
    for(i=0;i<q_field;i++){
        for(j=0;j<dc;j++)
            printf("%d ",result[i][j]);
        printf("\t");
        for(j=0;j<dc;j++)
            printf("%d ",R_Aux[i][j]);
        for(j=0;j<dc;j++){
            if(result[i][j] != R_Aux[i][j])
                error = 1;
        }
        printf("\n");
    }

    if(error ==1)
        printf("ERROR\n");
    else 
        printf("SUCCESS\n");

    asm volatile("jr x0;");
    return 0;

}
