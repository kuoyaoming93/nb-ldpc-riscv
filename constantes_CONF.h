/*
 * constantes.h
 *
 *  Created on: 22/09/2013
 *      Author: JesusOmar
 */

#ifndef CONSTANTES_H_
	#define CONSTANTES_H_

	#ifndef PI
		#define PI	3.14159265358979323846	/* pi */
	#endif

	/* Parametros del algoritmo */
	#define scaling_factor 0.9
	#define voto 1.5
	
	/* Parametros del decodificador */
    #define it_max 50 		/* Numero de iteraciones a simular para el algoritmo */
    #define Mchannels 1e10	/* Maximo numero de canales a simular */

   
    #define LLR_SCALE 1  /* 0 --> No Scale   ,   1 --> 1/sigma   ,    2 --> 2/sigma^2 */

   /* CUANTIFICACION PARA LOS BITS PROVENIENTES DEL CANAL*/
    #define RCH_QUANT 0
    #if RCH_QUANT == 1
        #define wRCH 4		/* Num bits totales para el RCH */
	    #define fRCH 3		/* Num bits fraccionales para el RCH */
	#endif

   /* CUANTIFICACION PARA LOS LLR */
    #define LLR_QUANT 0
    #if LLR_QUANT == 1
        #define wLN 8		/* Num bits totales para el LLR */
        #define fLN 4		/* Num bits fraccionales para el LLR */
	#endif

   /* CUANTIFICACION PARA EL DECODIFICADOR */
   #define DECO_QUANT 0
   #if DECO_QUANT == 1
	   #define wCN 7		/* Num bits totales para los mensajes de salida del check node */
	   #define fCN 3 		/* Num bits fraccionales para los mensajes de salida del check node */
	   #define wQN 8		/* Num bits totales para los mensajes entrantes al check node */
	   #define fQN 3 		/* Num bits fraccionales para los mensajes entrantes al check node */
   #else
      #define QN_SATp1 10000.0
	#endif

    
    /* Directorio para almacenar resultados */
	char resultspath[] =  "./RESULTS_ETMM_q16_dc4/";
	
#endif /* CONSTANTES_H_ */
