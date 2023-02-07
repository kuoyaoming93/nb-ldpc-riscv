/*
 ============================================================================
 Name        : main.c
 Author      : Jesus Lacruz, Yao-Ming Kuo
 Version     :
 Copyright   : 
 Description :
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "constantes_CONF_GF32.h"
#include "gf_tables.h"
#include "CNU_tables.h"
#include "Hmat_N837_M124_GF32.h"

//#define FOR_BACK

/*Functions*/
void awgn_channel(double r_ch[N_code][m_field], double sigma);
void LLR_function(double rx_points[N_code][m_field], double sigma, double constante, double Ln_aux[q_field][N_code]);
double randn_notrig(double mu, double sigma);

/* Constantes utilizadas para saturacion y truncado */
#if LLR_QUANT == 1
   double LLR_SATp1;
   double LLR_SATm1;
   double LLR_FLOORp1;
   double LLR_FLOORm1;
#endif
#if DECO_QUANT == 1
    double QN_SATp1;
    double QN_SATm1;
    double QN_FLOORp1;
    double QN_FLOORm1;

    double CN_SATp1;
    double CN_SATm1;
    double CN_FLOORp1;
    double CN_FLOORm1;
#endif

#if RCH_QUANT == 1
    double RCH_SATp1;
    double RCH_SATm1;
    double RCH_FLOORp1;
    double RCH_FLOORm1;
#endif

int wires[q_field][q_field] = {
    {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31},
		{1,0,19,6,30,11,3,28,23,21,17,5,20,24,15,14,25,10,31,2,12,9,26,8,13,16,22,29,7,27,4,18},
		{2,19,0,20,7,31,12,4,29,24,22,18,6,21,25,16,15,26,11,1,3,13,10,27,9,14,17,23,30,8,28,5},
		{3,6,20,0,21,8,1,13,5,30,25,23,19,7,22,26,17,16,27,12,2,4,14,11,28,10,15,18,24,31,9,29},
		{4,30,7,21,0,22,9,2,14,6,31,26,24,20,8,23,27,18,17,28,13,3,5,15,12,29,11,16,19,25,1,10},
		{5,11,31,8,22,0,23,10,3,15,7,1,27,25,21,9,24,28,19,18,29,14,4,6,16,13,30,12,17,20,26,2},
		{6,3,12,1,9,23,0,24,11,4,16,8,2,28,26,22,10,25,29,20,19,30,15,5,7,17,14,31,13,18,21,27},
		{7,28,4,13,2,10,24,0,25,12,5,17,9,3,29,27,23,11,26,30,21,20,31,16,6,8,18,15,1,14,19,22},
		{8,23,29,5,14,3,11,25,0,26,13,6,18,10,4,30,28,24,12,27,31,22,21,1,17,7,9,19,16,2,15,20},
		{9,21,24,30,6,15,4,12,26,0,27,14,7,19,11,5,31,29,25,13,28,1,23,22,2,18,8,10,20,17,3,16},
		{10,17,22,25,31,7,16,5,13,27,0,28,15,8,20,12,6,1,30,26,14,29,2,24,23,3,19,9,11,21,18,4},
		{11,5,18,23,26,1,8,17,6,14,28,0,29,16,9,21,13,7,2,31,27,15,30,3,25,24,4,20,10,12,22,19},
		{12,20,6,19,24,27,2,9,18,7,15,29,0,30,17,10,22,14,8,3,1,28,16,31,4,26,25,5,21,11,13,23},
		{13,24,21,7,20,25,28,3,10,19,8,16,30,0,31,18,11,23,15,9,4,2,29,17,1,5,27,26,6,22,12,14},
		{14,15,25,22,8,21,26,29,4,11,20,9,17,31,0,1,19,12,24,16,10,5,3,30,18,2,6,28,27,7,23,13},
		{15,14,16,26,23,9,22,27,30,5,12,21,10,18,1,0,2,20,13,25,17,11,6,4,31,19,3,7,29,28,8,24},
		{16,25,15,17,27,24,10,23,28,31,6,13,22,11,19,2,0,3,21,14,26,18,12,7,5,1,20,4,8,30,29,9},
		{17,10,26,16,18,28,25,11,24,29,1,7,14,23,12,20,3,0,4,22,15,27,19,13,8,6,2,21,5,9,31,30},
		{18,31,11,27,17,19,29,26,12,25,30,2,8,15,24,13,21,4,0,5,23,16,28,20,14,9,7,3,22,6,10,1},
		{19,2,1,12,28,18,20,30,27,13,26,31,3,9,16,25,14,22,5,0,6,24,17,29,21,15,10,8,4,23,7,11},
		{20,12,3,2,13,29,19,21,31,28,14,27,1,4,10,17,26,15,23,6,0,7,25,18,30,22,16,11,9,5,24,8},
		{21,9,13,4,3,14,30,20,22,1,29,15,28,2,5,11,18,27,16,24,7,0,8,26,19,31,23,17,12,10,6,25},
		{22,26,10,14,5,4,15,31,21,23,2,30,16,29,3,6,12,19,28,17,25,8,0,9,27,20,1,24,18,13,11,7},
		{23,8,27,11,15,6,5,16,1,22,24,3,31,17,30,4,7,13,20,29,18,26,9,0,10,28,21,2,25,19,14,12},
		{24,13,9,28,12,16,7,6,17,2,23,25,4,1,18,31,5,8,14,21,30,19,27,10,0,11,29,22,3,26,20,15},
		{25,16,14,10,29,13,17,8,7,18,3,24,26,5,2,19,1,6,9,15,22,31,20,28,11,0,12,30,23,4,27,21},
		{26,22,17,15,11,30,14,18,9,8,19,4,25,27,6,3,20,2,7,10,16,23,1,21,29,12,0,13,31,24,5,28},
		{27,29,23,18,16,12,31,15,19,10,9,20,5,26,28,7,4,21,3,8,11,17,24,2,22,30,13,0,14,1,25,6},
		{28,7,30,24,19,17,13,1,16,20,11,10,21,6,27,29,8,5,22,4,9,12,18,25,3,23,31,14,0,15,2,26},
		{29,27,8,31,25,20,18,14,2,17,21,12,11,22,7,28,30,9,6,23,5,10,13,19,26,4,24,1,15,0,16,3},
		{30,4,28,9,1,26,21,19,15,3,18,22,13,12,23,8,29,31,10,7,24,6,11,14,20,27,5,25,2,16,0,17},
		{31,18,5,29,10,2,27,22,20,16,4,19,23,14,13,24,9,30,1,11,8,25,7,12,15,21,28,6,26,3,17,0}
};

int main(int argc, char * argv[]) {

	int i,j,k,eS,eC,var,row; /* Used for loop indexs */

	double SNRdB[10];
	double sigma[10];
	double sigma_2[10];

	int 	decoded[N_code][m_field];
	double 	r_ch[N_code][m_field];
	double 	Ln_aux[q_field][N_code];
	double 	Rmn_SRL[M_code][q_field][dc];

	double Qmn[q_field][dc];
	double Qmn_temp[q_field];

	int aux1, aux2;
	double MAX_temp; //MAX_temp2;

    double dQ_min1, dQ_min2; 
    int dQ_pos1, dQ_pos2;
    
	int beta;
	int z[dc];
	double dWmn[q_field][dc];
	int temp[q_field][dc];
	double min1[q_field], min2[q_field];
	int pos[q_field];
	double max[q_field/2];
	int cam1_temp[q_field/2 - 1], cam2_temp[q_field/2 - 1];
	double min_global[q_field]; //min2_global[q_field];
	int cam1[q_field], cam2[q_field];
	double Rmn[q_field][dc];
	double Qn_NEW[q_field][dc];

#ifdef FOR_BACK
	int a;
	double R_Forward[q_field][dc], R_Backward[q_field][dc], R_Aux[q_field][dc];
    double R_aux[q_field], R_aux_2[q_field], R_compare[q_field], R_compare_2[q_field], R_Backward_aux[q_field];
    double min_temp, min_temp_2;
#endif

	int H_decoded[it_max];
	int MNBE_Hdecoded[it_max];
	int MNPE_Hdecoded[it_max];

	int error = 0;

	char nombre_file[120];
	//char nombre_file_rng[120];
	FILE* archivo;

	int LPerrors = 0;		/* Numero de errores a buscar */
	int semilla = 0;
	double EbNodB[10];
	int EbNo_NUM = 0;
	
   #if LLR_QUANT == 1
       LLR_SATp1 = pow(2,wLN-fLN-1)-pow(2,-fLN);
	    LLR_SATm1 = -pow(2,wLN-fLN-1);
	    LLR_FLOORp1 = pow(2,fLN);
	    LLR_FLOORm1 = pow(2,-fLN);
   #endif

	 #if DECO_QUANT == 1
	    QN_SATp1 = pow(2,wQN-fQN-1)-pow(2,-fQN);
	    QN_SATm1 = -pow(2,wQN-fQN-1);
	    QN_FLOORp1 = pow(2,fQN);
	    QN_FLOORm1 = pow(2,-fQN);

	    CN_SATp1 = pow(2,wCN-fCN-1)-pow(2,-fCN);
	    CN_SATm1 = -pow(2,wCN-fCN-1);
	    CN_FLOORp1 = pow(2,fCN);
	    CN_FLOORm1 = pow(2,-fCN);

    #endif

    #if RCH_QUANT == 1
	    RCH_SATp1 = pow(2,wRCH-fRCH-1)-pow(2,-fRCH);
	    RCH_SATm1 = -pow(2,wRCH-fRCH-1);
	    RCH_FLOORp1 = pow(2,fRCH);
	    RCH_FLOORm1 = pow(2,-fRCH);
    #endif

	if (argc < 4)
	{
      puts("No se han pasado suficientes argumentos al programa...");
		return EXIT_SUCCESS;
	}
  	else
	{
	    semilla = atoi(argv[1]);
        LPerrors = atoi(argv[2]);
		EbNo_NUM = argc - 3;
        for (i=0; i < EbNo_NUM; i++)
		{
			EbNodB[i] = atof(argv[i+3]);
		}
	}

	for (i=0; i<EbNo_NUM;i++)
	{
		SNRdB[i] = 10*log10(2*rate*pow(10,EbNodB[i]/10));
		sigma[i] = sqrt(pow(10,(-SNRdB[i]/10)));
      
        #if LLR_SCALE == 2		
            sigma_2[i] = (1/(sigma[i]*sigma[i]*2));
        #endif		
        #if LLR_SCALE == 1
            sigma_2[i] = 1/sigma[i];
        #endif
        #if LLR_SCALE == 0
            sigma_2[i] = 1;  
        #endif		
	}

	/* srand(time(0)); */
	srand(semilla);

            /* Loop through all EbNodB values */
            for (eS=0; eS<EbNo_NUM;eS++)
            {

                sprintf(nombre_file, "%s(%u,%u)_L%u[%u]_N%u_K%u_EBN%4.2f.txt", resultspath, N_code, K_code, LPerrors, semilla, Nbpb, (int) (Nbpb*rate), EbNodB[eS]);
		        //sprintf(nombre_file_rng, "%s(%u,%u)_L%u[%u]_N%u_K%u_EBN%4.2f_RNG", resultspath, N_code, K_code, LPerrors, semilla, Nbpb, (int) (Nbpb*rate), EbNodB[eS]);  
				
		        error = 0;
		        archivo = fopen(nombre_file, "r");
		        if (archivo!= NULL)
		        {
	                char s1[] = "MNPE_Hdecoded";
			        char ss1[15];
			        char s2[] = "MNBE_Hdecoded";
			        char s3[] = "eC";

			        /*while (!feof(archivo))
			        {*/
			        if (fscanf(archivo,"%s",ss1) == 1)
			        {
				        if (strcmp(s1,ss1)==0)
				        {

					        for (i=0; i<it_max; i++)
					        {
						        if (fscanf(archivo,"%d",&MNPE_Hdecoded[i]) != 1)
							        error = 1;
					        }
				        }
				        else
					        error = 1;
			        }
			        else
				        error = 1;

			        if (fscanf(archivo,"%s",ss1) == 1)
			        {
				        if (strcmp(s2,ss1)==0)
				        {
					        for (i=0; i<it_max; i++)
					        {
						        if (fscanf(archivo,"%d",&MNBE_Hdecoded[i]) != 1)
							        error = 1;
					        }
				        }
				        else
					        error = 1;

			        }
			        else
				        error = 1;

		            if (fscanf(archivo,"%s",ss1) == 1)
			        {
				        if (strcmp(s3,ss1)==0)
				        {
					        if (fscanf(archivo,"%d",&eC) != 1)
						        error = 1;
				        }
				        else
					        error = 1;
			        }
			        else
				        error = 1;
			        /*}*/

			        fclose(archivo);

			        if ((MNPE_Hdecoded[it_max-1] == LPerrors) && (error == 0))
				        continue;


			        /* Lee el estado del generador aleatorio */
			        /*
			        FILE * file = fopen(nombre_file_rng, "r");

			        if (!file)
				        error = 1;
			        else
			        {
				        gsl_rng_fread(file, rng);
				        fclose(file);
			        }
                    */
			        /***************************************/

		        }
		        else
			        error = 1;

		        if (error == 1)
		        {
			        eC = 0;
			        for (i=0; i<it_max; i++)
			        {
				        MNBE_Hdecoded[i] =0;
				        MNPE_Hdecoded[i] = 0;
			        }
			        //gsl_rng_set(rng, semilla);
		        }

		

		        while (eC < Mchannels)
		        {
			        /* Passing signal through the AWGN channel */
			        awgn_channel(r_ch, sigma[eS]);

			        /* LLRs_function Ln(a)=ln(Pr(cn=zn|channel)/Pr(cn=a|channel)) */
			        LLR_function(r_ch, sigma[eS], sigma_2[eS], Ln_aux);

					//if(eC == 1){
						/*printf("Message with noise: ");
						for (i = 0; i< N_code; i++)
							for(j=0; j<m_field;j++)
								printf("r_ch[%d][%d]: %f\n",i,j,r_ch[i][j]);*/

						/*printf("LLR Function: ");
						for (i = 0; i< q_field; i++){
							for(j=0; j<N_code;j++)
								printf("%f, \n",Ln_aux[i][j]);
							printf("\n");
						}*/
//								printf("Ln_aux[%d][%d]: %f\n",i,j,Ln_aux[i][j]);
					//}
					
					
			        /* Inicialize Rmn_SRL to zero values */
			        for (i = 0; i<M_code; i++)
				        for (j=0; j<q_field; j++)
					        for (k=0; k<dc; k++)
						        Rmn_SRL[i][j][k] = 0.0;

			        for (var=0; var<it_max; var++)
				        H_decoded[var] = 0;

			         /* Loop through all decoding iterations */
			         for (var=0; var<it_max; var++)
			         {						
				        /* Loop through all rows of parity check matrix */
				        for (row=0; row<M_code; row++)
                        {

                            beta = 0; /* Inicializa el valor del sindrome a cero */
                            /* Extract Qmn(a) messages from Qn(a) memories */
				            for (i=0;i<dc;i++)
	                        {

					            for (j=0;j<q_field;j++)
                                {
                                    Qmn_temp[j] = Ln_aux[j][col[row][i]];
	                            }

						        /* PERMUTACION DE LA COLUMNA EXTRAIDA */

						        if (pow_coefH[row][i]==0)
		                        {
		                            for (j=0;j<q_field;j++)
		                            {
		                                Qmn[j][i] = Qmn_temp[j];
		                            }
		                        }
		                        else
		                        {
		                            aux2 = pow_coefH[row][i] + 1;
		                            aux1 = q_field - aux2;

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
                                        z[i] = pos_exp_table[j];  
                                    }
                                }
                                beta = gfadd[beta][z[i]];  /* SINDROME DEL CHECK NODE */

	                            /*****************************************/

	                            /* Normalizacion y precision finita*/
	                            for (j=0;j<q_field;j++)
	                            {
						            #if DECO_QUANT == 1
                                    Qmn[j][i] = floor((Qmn[j][i]-MAX_temp)*QN_FLOORp1)*QN_FLOORm1; 
					                    if(Qmn[j][i] > QN_SATp1)
						                    Qmn[j][i] = QN_SATp1;
							        #else
                                        Qmn[j][i] = Qmn[j][i]-MAX_temp; 
					                #endif
	                            }
                                /*****************/
	                        }

                            /*******************************************/
                            /**********  CNU FUNCTION ******************/
                            /*******************************************/
#ifdef FOR_BACK
							// Clean Forward and Backward variables and set value to first and last column
							// Clean R_Aux
							for(j=0;j<dc;j++){
								for(i=0;i<q_field;i++){        
									if(j==0)
										R_Forward[i][j] = Qmn[i][j];
									else
										R_Forward[i][j] = 0;

									if(j==dc-1){
										R_Backward[i][j] = Qmn[i][j];
										R_Aux[i][j] = Qmn[i][j];
									}
									else{
										R_Backward[i][j] = 0;
										R_Aux[i][j] = 0;
									}            
								}
							}

							for(a=2;a<dc;a++){
								for(i=0;i<q_field;i++){
									R_aux[i]    = Qmn[i][a-1];
									R_aux_2[i]  = Qmn[i][4-a];
								}
								for(i=0;i<q_field;i++){

									// Search maximum
									for(k=0;k<q_field;k++){

										// Forward
										if(R_Forward[k][a-2] > R_aux[wires[k][i]])
											R_compare[k] = R_Forward[k][a-2];
										else
											R_compare[k] = R_aux[wires[k][i]];

										// Backward
										if(R_Backward[k][5-a] > R_aux_2[wires[k][i]])
											R_compare_2[k] = R_Backward[k][5-a];
										else
											R_compare_2[k] = R_aux_2[wires[k][i]];
									}

									// Search minimum
									min_temp = 10000;
									min_temp_2 = 10000;
									for(k=0;k<q_field;k++){
										// Forward
										if(R_compare[k] < min_temp)
											min_temp = R_compare[k];
										// Backward
										if(R_compare_2[k] < min_temp_2)
											min_temp_2 = R_compare_2[k];
									}

									R_Forward[i][a-1] = min_temp;
									R_Backward[i][4-a] = min_temp_2;
								}
							}

							// Standard Min Max
							for(a=2;a<dc;a++){
								for(i=0;i<q_field;i++)
									R_Backward_aux[i] = R_Backward[i][a];

								for(i=0;i<q_field;i++){
									// Search maximum
									for(k=0;k<q_field;k++){
										if(R_Forward[k][a-2] > R_Backward_aux[wires[k][i]])
												R_compare[k] = R_Forward[k][a-2];
											else
												R_compare[k] = R_Backward_aux[wires[k][i]];
									}
									// Search minimum
									min_temp = 10000;
									for(k=0;k<q_field;k++){
										if(R_compare[k] < min_temp)
											min_temp = R_compare[k];
									}

									Rmn[i][a-1] = min_temp;
								}        
							}

							for(i=0;i<q_field;i++){
								Rmn[i][0] =       R_Backward[i][1];
								Rmn[i][dc-1] =    R_Forward[i][dc-2];
							}
#else                            
				            /* paso de mensajes al dominio delta */
					        for (i=0;i<dc;i++)
					        {
					            for (j=0;j<q_field;j++)
					            {
					                temp[j][i] = gfadd[z[i]][pos_exp_table[j]];
							        dWmn[j][i] = Qmn[ exp_pos_table[ temp[j][i] ] -1][i];
						        }
					        }


					        /* Busqueda de minimos */
					        for (j=0;j<q_field;j++)
					        {
			                    double min_temp = dWmn[j][0];
						        double min_temp2 = 10000;
						        int pos_temp = 0;
						   
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
	                            double min_temp = min1[i+1];
                                int pos_temp = 0;
                                
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
                            
                            /* TEMPORAL PRINTF */
                            /*
                            printf("dQ = ");
                            for (i = 0; i< q_field; i++)
                                printf("%7.2f ",min_global[i]);
                            printf("\n");*/
                            
                            
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

                            /*
                            printf("dQ_min1 = %7.2f       pos = %d \n", dQ_min1,dQ_pos1);
                            printf("dQ_min2 = %7.2f       pos = %d \n", dQ_min2,dQ_pos2);
                            return;
                            */
#endif						   
                            
                            for (i=0;i<dc;i++)
					        {
#ifndef FOR_BACK								
				                double oRmn_temp;
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
								            oRmn_temp = voto * dQ_min2;
							        }
							        
						            temp[j][i] = exp_pos_table[gfadd[temp[j][i]][beta]]-1;
						         
                                    /*Aplica escalado y precision finita a los mensajes de salida del CN*/
                                    #if DECO_QUANT == 1						   
                                        Rmn[temp[j][i]][i] = floor((oRmn_temp*scaling_factor)*CN_FLOORp1)*CN_FLOORm1; 
			                                if(Rmn[temp[j][i]][i] > CN_SATp1)
					                            Rmn[temp[j][i]][i] = CN_SATp1;
				                    #else
                                        Rmn[temp[j][i]][i] = oRmn_temp*scaling_factor;
                                    #endif
					            }
#endif
					            for (j=0;j<q_field;j++)
					            {
						        
                                    Rmn_SRL[row][j][i] = Rmn[j][i];

						            /*CNout = CNout + Qmn_prima;*/
						            Rmn[j][i] = Rmn[j][i] + Qmn[j][i];
					            }

					            /****** PERMUTACION INVERSA *******/
						        aux1 = pow_coefH[row][i];
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
						        /***********************************/

					            /* PRECISION FINITA Y ACTUALIZACION DE LAS VN*/

					            for (j=0;j<q_field;j++)
					            {
                                    #if DECO_QUANT == 1
					                    Ln_aux[j][col[row][i]] = floor((Qn_NEW[j][i])*QN_FLOORp1)*QN_FLOORm1; 
						                if(Ln_aux[j][col[row][i]] > QN_SATp1)
						                    Ln_aux[j][col[row][i]] = QN_SATp1;
				                    #else
                                        Ln_aux[j][col[row][i]] = Qn_NEW[j][i];
                                    #endif
					            }

				            }

					        /******* TEMPORAL *************/
					        
					        /*for (i=0; i<dc; i++)
					        {
						        for (j=0; j<q_field; j++)
						        {
							        printf("%7.4f ", Qn_NEW[j][i]);
						        }
						        printf("\n");
					        }
					        printf("\n");*/
					        
					        /******************************/

				        }

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

				        H_decoded[var] = aux2;

				        if (aux2 == 0)
					        break;

			        }

			        /********* TEMPORAL **************/
			        
			        /*printf("\n Hdecoded = ");
			        for (i=0; i<it_max; i++)
			        {
				        printf("%d ",H_decoded[i]);
			        }*/
			        

			        for (i=0; i<it_max; i++)
			        {
				        MNBE_Hdecoded[i] = MNBE_Hdecoded[i] + H_decoded[i];
				        MNPE_Hdecoded[i] = MNPE_Hdecoded[i] + (int)(H_decoded[i]>0);
			        }

			        if (MNPE_Hdecoded[it_max-1] == LPerrors)
				        break;

			        eC++;

			        if ((eC%10000)==1)
			        {
				        /*
				        FILE * file = fopen(nombre_file_rng, "w");
				        gsl_rng_fwrite(file, rng);
				        fclose(file);
                        */
				        archivo = fopen(nombre_file, "w");
				        if (archivo)
				        {
					        fprintf(archivo,"MNPE_Hdecoded\n");
					        for (i=0; i<it_max; i++)
					        {
						        fprintf(archivo,"%d \n",MNPE_Hdecoded[i]);
					        }
					        fprintf(archivo,"MNBE_Hdecoded\n");
					        for (i=0; i<it_max; i++)
					        {
						        fprintf(archivo,"%d\n",MNBE_Hdecoded[i]);
					        }
					        fprintf(archivo,"eC\n%d \n",eC);

					        fclose(archivo);
				        }			
				        
				        printf("\n MNPE_Hdecoded = ");
				        for (i=0; i<it_max; i++)
				        {
					        printf("%d ", MNPE_Hdecoded[i]);
				        }
				        printf("\n MNBE_Hdecoded = ");
				        for (i=0; i<it_max; i++)
				        {
					        printf("%d ", MNBE_Hdecoded[i]);
				        }
				        printf("\n eC = %d \n",eC);
				        
			        }
		        }
                /*
		        FILE * file = fopen(nombre_file_rng, "w");
		        gsl_rng_fwrite(file, rng);
		        fclose(file);
                */
		        archivo = fopen(nombre_file, "w");
		        if (archivo)
		        {
			        fprintf(archivo,"MNPE_Hdecoded\n");
			        for (i=0; i<it_max; i++)
			        {
				        fprintf(archivo,"%d \n",MNPE_Hdecoded[i]);
			        }
			        fprintf(archivo,"MNBE_Hdecoded\n");
			        for (i=0; i<it_max; i++)
			        {
				        fprintf(archivo,"%d\n",MNBE_Hdecoded[i]);
			        }
			        fprintf(archivo,"eC\n%d \n",eC);

			        fclose(archivo);
		        }
		        /*
		        remove(nombre_file_rng);
                */
    }

	return EXIT_SUCCESS;
}



void awgn_channel(double r_ch[N_code][m_field], double sigma)
{
	int i,j;

	for (i = 0; i< N_code; i++)
	{
		for (j=0; j<m_field; j++)
		{
			#if RCH_QUANT == 1
        		//r_ch[i][j] = floor( ((double) m[i][j] + gsl_ran_gaussian(rng,sigma)) * RCH_FLOORp1) * RCH_FLOORm1; 
				r_ch[i][j] = floor( ((double) m[i][j] + randn_notrig(0.0,sigma)) * RCH_FLOORp1) * RCH_FLOORm1; 
				if(r_ch[i][j] > RCH_SATp1)
					r_ch[i][j] = RCH_SATp1;
				if (r_ch[i][j] < RCH_SATm1) 
					r_ch[i][j] = RCH_SATm1;
            #else
                //r_ch[i][j] = (double) m[i][j] + gsl_ran_gaussian(rng,sigma);
                r_ch[i][j] = (double) m[i][j] + randn_notrig(0.0,sigma);
            #endif
		}
	}
}

void LLR_function(double rx_points[N_code][m_field], double sigma, double constante, double Ln_aux[q_field][N_code])
{
   int i,j,k;
   double MIN_temp;

    for (i=0;i<N_code;i++)
    {
        MIN_temp = 10000;
        for (j=0;j<q_field;j++)
        {
            Ln_aux[j][i]=0.0;

            for (k=0;k<m_field;k++)
            {
                Ln_aux[j][i]=Ln_aux[j][i]+constante*rx_points[i][k]*fld_bpsk_order[j][k];
            }

			/* Busca el minimo de los LLR para normalizar los mismos */
            if (MIN_temp > Ln_aux[j][i])
            {
                MIN_temp = Ln_aux[j][i];
            }
        }

        /* Normaliza los LLR calculados y los cuantifica*/
        for (j=0;j<q_field;j++)
        {
            #if LLR_QUANT == 1
           		Ln_aux[j][i] = floor((Ln_aux[j][i]-MIN_temp)*LLR_FLOORp1)*LLR_FLOORm1; 
				   if(Ln_aux[j][i] > LLR_SATp1)
					   Ln_aux[j][i] = LLR_SATp1;
				   if (Ln_aux[j][i]< LLR_SATm1) 
					   Ln_aux[j][i] = LLR_SATm1;
            #else
               Ln_aux[j][i] = (Ln_aux[j][i]-MIN_temp);
            #endif
        }
    }

    /* TEMP PRINTF */
    /*   
    printf("LLR = \n");
    for (i=0;i<N_code;i++)
    {
        for(j=0;j<q_field;j++)
        {
            printf("%6.4f, ",Ln_aux[j][i]);
        }
        printf("\n");
    }
    */

    return;
}

/*
void save_state(char nombre[100])
{
	FILE * file = fopen(nombre, "w");
	gsl_rng_fwrite(file, rng);
	fclose(file);
}
*/

//Channel
//double randn_notrig(double mu=0.0, double sigma=1.0) //http://www.dreamincode.net/code/snippet1446.htm
double randn_notrig(double mu, double sigma) //http://www.dreamincode.net/code/snippet1446.htm
{
    static bool deviateAvailable = false; // flag
    static float storedDeviate;   // deviate from previous calculation
    double polar, rsquared, var1, var2;
 
 
    // If no deviate has been stored, the polar Box-Muller transformation is 
    // performed, producing two independent normally-distributed random
    // deviates.  One is stored for the next round, and one is returned.
    if (!deviateAvailable) {

        // choose pairs of uniformly distributed deviates, discarding those 
        // that don't fall within the unit circle
        do 
        {
            //var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
            //var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
            var1 = 2.0 * ( ((double) rand()) / ((double) RAND_MAX) ) - 1.0;
            var2 = 2.0 * ( ((double) rand()) / ((double) RAND_MAX) ) - 1.0;
            rsquared=var1*var1+var2*var2;
        } 
        while ( rsquared>=1.0 || rsquared == 0.0);

        // calculate polar tranformation for each deviate
        polar=sqrt(-2.0*log(rsquared)/rsquared);
  
        // store first deviate and set flag
        storedDeviate=var1*polar;
        deviateAvailable=true;
  
        // return second deviate
        return var2*polar*sigma + mu;
    }
 
    // If a deviate is available from a previous call to this function, it is
    // returned, and the flag is set to false.
    else 
    {
        deviateAvailable=false;
        return storedDeviate*sigma + mu;
    }
}
//end Channel

