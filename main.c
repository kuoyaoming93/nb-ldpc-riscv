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

#include "constantes_CONF.h"
#include "gf_tables_GF16.h"
#include "CNU_tables_GF16.h"
#include "Hmat_N32_M16_GF16.h"


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
            sigma_2[i] = (2/(sigma[i]*sigma[i]));
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
                            
                            for (i=0;i<dc;i++)
					        {
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

