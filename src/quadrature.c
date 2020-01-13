/*
 * quadrature.c
 *
 *  Created on: July 2, 2014
 *      Author: yukikoshimizu
 */

#include "quadrature.h"

void quadrature(GuassianQuad quadra){
	int n=quadra.n;
	if (n==1){
		quadra.w.array[0] = 2; quadra.x_i.array[0]=0;
	}
	else if (n==2){
		quadra.w.array[0]=1.0;	quadra.x_i.array[0]=-0.5773502691896257;
		quadra.w.array[1]=1.0;	quadra.x_i.array[1]=0.5773502691896257;
	}
	else if (n==3){
		quadra.w.array[0]=0.8888888888888888;	quadra.x_i.array[0]=0;
		quadra.w.array[1]=0.5555555555555554;	quadra.x_i.array[1]=-0.7745966692414834;
		quadra.w.array[2]=0.5555555555555554;	quadra.x_i.array[2]=0.7745966692414834;
	}
	else if (n==4){
		quadra.w.array[0]=0.6521451548625461;	quadra.x_i.array[0]=-0.3399810435848563;
		quadra.w.array[1]=0.6521451548625461;	quadra.x_i.array[1]=0.3399810435848563;
		quadra.w.array[2]=0.3478548451374538;	quadra.x_i.array[2]=-0.8611363115940526;
		quadra.w.array[3]=0.3478548451374538;	quadra.x_i.array[3]=0.8611363115940526;
	}
	else if (n==5){
		quadra.w.array[0]=0.5688888888888889;	quadra.x_i.array[0]=0;
		quadra.w.array[1]=0.4786286704993665;	quadra.x_i.array[1]=-0.5384693101056831;
		quadra.w.array[2]=0.4786286704993665;	quadra.x_i.array[2]=0.5384693101056831;
		quadra.w.array[3]=0.2369268850561891;	quadra.x_i.array[3]=-0.9061798459386640;
		quadra.w.array[4]=0.2369268850561891;	quadra.x_i.array[4]=0.9061798459386640;
	}
	else if (n==6){
		quadra.w.array[0]=0.3607615730481386;	quadra.x_i.array[0]=0.6612093864662645;
		quadra.w.array[1]=0.3607615730481386;	quadra.x_i.array[1]=-0.6612093864662645;
		quadra.w.array[2]=0.4679139345726910;	quadra.x_i.array[2]=-0.2386191860831969;
		quadra.w.array[3]=0.4679139345726910;	quadra.x_i.array[3]=0.2386191860831969;
		quadra.w.array[4]=0.1713244923791704;	quadra.x_i.array[4]=-0.9324695142031521;
		quadra.w.array[5]=0.1713244923791704;	quadra.x_i.array[5]=0.9324695142031521;
	}
	else if (n==7){
		quadra.w.array[0]=0.4179591836734694;	quadra.x_i.array[0]=0;
		quadra.w.array[1]=0.3818300505051189;	quadra.x_i.array[1]=0.4058451513773972;
		quadra.w.array[2]=0.3818300505051189;	quadra.x_i.array[2]=-0.4058451513773972;
		quadra.w.array[3]=0.2797053914892766;	quadra.x_i.array[3]=-0.7415311855993945;
		quadra.w.array[4]=0.2797053914892766;	quadra.x_i.array[4]=0.7415311855993945;
		quadra.w.array[5]=0.1294849661688697;	quadra.x_i.array[5]=-0.9491079123427585;
		quadra.w.array[6]=0.1294849661688697;	quadra.x_i.array[6]=0.9491079123427585;
	}
	else if (n==8){
		quadra.w.array[0]=0.3626837833783620;	quadra.x_i.array[0]=-0.1834346424956498;
		quadra.w.array[1]=0.3626837833783620;	quadra.x_i.array[1]=0.1834346424956498;
		quadra.w.array[2]=0.3137066458778873;	quadra.x_i.array[2]=-0.5255324099163290;
		quadra.w.array[3]=0.3137066458778873;	quadra.x_i.array[3]=0.5255324099163290;
		quadra.w.array[4]=0.2223810344533745;	quadra.x_i.array[4]=-0.7966664774136267;
		quadra.w.array[5]=0.2223810344533745;	quadra.x_i.array[5]=0.7966664774136267;
		quadra.w.array[6]=0.1012285362903763;	quadra.x_i.array[6]=-0.9602898564975363;
		quadra.w.array[7]=0.1012285362903763;	quadra.x_i.array[7]=0.9602898564975363;
	}
	else if (n==9){
		quadra.w.array[0]=0.3302393550012598;	quadra.x_i.array[0]=0;
		quadra.w.array[1]=0.1806481606948574;	quadra.x_i.array[1]=-0.8360311073266358;
		quadra.w.array[2]=0.1806481606948574;	quadra.x_i.array[2]=0.8360311073266358;
		quadra.w.array[3]=0.0812743883615744;	quadra.x_i.array[3]=-0.9681602395076261;
		quadra.w.array[4]=0.0812743883615744;	quadra.x_i.array[4]=0.9681602395076261;
		quadra.w.array[5]=0.3123470770400029;	quadra.x_i.array[5]=-0.3242534234038089;
		quadra.w.array[6]=0.3123470770400029;	quadra.x_i.array[6]=0.3242534234038089;
		quadra.w.array[7]=0.2606106964029354;	quadra.x_i.array[7]=-0.6133714327005904;
		quadra.w.array[8]=0.2606106964029354;	quadra.x_i.array[8]=0.6133714327005904;
	}
	else if (n==10){
		quadra.w.array[0]=0.2955242247147529;	quadra.x_i.array[0]=-0.1488743389816312;
		quadra.w.array[1]=0.2955242247147529;	quadra.x_i.array[1]=0.1488743389816312;
		quadra.w.array[2]=0.2692667193099963;	quadra.x_i.array[2]=-0.4333953941292472;
		quadra.w.array[3]=0.2692667193099963;	quadra.x_i.array[3]=0.4333953941292472;
		quadra.w.array[4]=0.2190863625159820;	quadra.x_i.array[4]=-0.6794095682990244;
		quadra.w.array[5]=0.2190863625159820;	quadra.x_i.array[5]=0.6794095682990244;
		quadra.w.array[6]=0.1494513491505806;	quadra.x_i.array[6]=-0.8650633666889845;
		quadra.w.array[7]=0.1494513491505806;	quadra.x_i.array[7]=0.8650633666889845;
		quadra.w.array[8]=0.0666713443086881;	quadra.x_i.array[8]=-0.9739065285171717;
		quadra.w.array[9]=0.0666713443086881;	quadra.x_i.array[9]=0.9739065285171717;
	}
}
