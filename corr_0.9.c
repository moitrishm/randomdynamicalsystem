#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_sf_hyperg.h>


int main (void)
{
	double i, prob;
	float m, temp, x_initial;
	char stn[40];
	char stn1[40];
	char stnp[40];

	float results[50] = { 0.0000000000 };
	int numbi = 100000;
	float xe[1] = { 0.0000000000 };
	time_t t0, t00, t0p;
	long t1, t11, t1p;
	double p=0.9;
	int interval=0;
	double probs[20] = { 0.0000000000 };
	double values[20] = { 0.0000000000 };


	mpf_t x, probm;
	mp_bitcnt_t n=256;
	mpf_init2(x,n);
	mpf_init2(probm,n);
	mpz_t seed, seed1, seedp;
	
	
	gmp_randstate_t rstate, rstate1, rstatep;
	gmp_randinit_mt (rstate);
	gmp_randinit_mt (rstate1);
	gmp_randinit_mt (rstatep);
	mpfr_t s, t, t_init, start, init, t_final;
	mpfr_init2 (t, 1000);

	mpfr_init2 (s, 1000);
	mpfr_init2 (init, 1000);

	mpfr_init2 (start, 1000);
	mpfr_init2 (t_init, 1000);
	mpfr_init2 (t_final, 1000);

	mpfr_set_d(t_init, 0, MPFR_RNDD);
	mpfr_set_d(t_final, 0, MPFR_RNDD);

	//xs=0.18367;
	//ys=0.082997;	

  //mpfr_set_d (s, 1.0, MPFR_RNDD);
	//mpfr_init2 (u, 1000);
	//mpfr_set_d (u, k, MPFR_RNDD);
	FILE *gnuplot_9u_cor = fopen("gnuplot_9u_cor.txt", "w");

	for( ;interval <20; interval++){
		probs[interval] = pow(0.5, (interval+1));
	}

	values[0] = ((2*p - 1)/(3*p - 2)) * (1 - pow((2*(1-p)/p) , 1)) * probs[0];

	for(interval = 1; interval < 20; interval++){
		values[interval] =  values[interval-1] + (((2*p - 1)/(3*p - 2)) * (1 - pow((2*(1-p)/p) , (interval+1)))*probs[interval] );

	}
	
	

  //fprintf(gnuplot, "plot '-'\n");
  
	for(int l=0; l<numbi; l++)
	{
		t0=time(NULL);
		t1=t0*(l+1);
		sprintf(stn, "%ld", t1);
		mpz_init_set_str(seed, stn,10);
		gmp_randseed (rstate, seed);
		mpfr_urandomb(init,rstate);
		t0p=time(NULL);
		t1p=t0p*(l+1);
		sprintf(stnp, "%ld", t1p);
		mpz_init_set_str(seedp, stnp,10);			
		gmp_randseed (rstatep, seedp);
		mpf_urandomb(probm,rstatep,n);
		prob = mpf_get_d(probm);
		//printf("%f\n", prob);
		if (prob <= values[0]){

			mpfr_set_d(start, probs[0], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[0] && prob <= values[1]){

			mpfr_set_d(start, probs[1], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);

		}
		else if (prob > values[1] && prob <= values[2]){

			mpfr_set_d(start, probs[2], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob > values[2] && prob <= values[3]){

			mpfr_set_d(start, probs[3], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}

		else if (prob >values[3] && prob <= values[4]){
			
			mpfr_set_d(start, probs[4], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[4] && prob <= values[5]){
			
			mpfr_set_d(start, probs[5], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[5] && prob <= values[6]){
			
			mpfr_set_d(start, probs[6], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[6] && prob <= values[7]){
			
			mpfr_set_d(start, probs[7], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[7] && prob <= values[8]){
			
			mpfr_set_d(start, probs[8], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[8] && prob <= values[9]){
			
			mpfr_set_d(start, probs[9], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[9] && prob <= values[10]){
			
			mpfr_set_d(start, probs[10], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[10] && prob <= values[11]){
			
			mpfr_set_d(start, probs[11], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[11] && prob <= values[12]){
			
			mpfr_set_d(start, probs[12], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[12] && prob <= values[13]){
			
			mpfr_set_d(start, probs[13], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[13] && prob <= values[14]){
			
			mpfr_set_d(start, probs[14], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[14] && prob <= values[15]){
			
			mpfr_set_d(start, probs[15], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[15] && prob <= values[16]){
			
			mpfr_set_d(start, probs[16], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[16] && prob <= values[17]){
			
			mpfr_set_d(start, probs[17], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[17] && prob <= values[18]){
			
			mpfr_set_d(start, probs[18], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}
		else if (prob >values[18] && prob <= values[19]){
			
			mpfr_set_d(start, probs[19], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add(t_init, start, init, MPFR_RNDD);
		}

		else if(prob > values[19]){
			mpfr_set_d(start, probs[19], MPFR_RNDD);
			mpfr_mul(init, init, start, MPFR_RNDD);
			mpfr_add_d(t_init, init, 0, MPFR_RNDD);

		}
	

		
		//mpfr_out_str (stdout, 10, 10, t_init, MPFR_RNDD);
		//putchar ('\n');
		mpfr_set (t, t_init, MPFR_RNDD);
		x_initial=mpfr_get_flt(t_init, MPFR_RNDD);
		xe[0] = xe[0] + x_initial;
		//printf("%.10f\n", x_initial );
		mpfr_mul(t_final, t, t, MPFR_RNDD);
		m= mpfr_get_flt (t_final, MPFR_RNDD);
		//printf("%.10f\n", m);
		results[0] = results[0] + m;
		for (int j = 1; j < 50; j++)
		{

			//mpfr_set(t_init, t, MPFR_RNDD);
			t00=time(NULL);
			t11=t00*(j+1);
			sprintf(stn1, "%ld", t11);
			mpz_init_set_str(seed1, stn1,10);			
			gmp_randseed (rstate1, seed1);
			mpf_urandomb(x,rstate1,n);
			
			i = mpf_get_d(x);
			//mpfr_set (t, u, MPFR_RNDD);
			if ( i <= p){
				mpfr_set_d (s, 2.0, MPFR_RNDD);	
				mpfr_mul(t, t, s, MPFR_RNDU);
				mpfr_frac (t, t, MPFR_RNDD);
 
			}
			else{
				mpfr_set_d (s, 0.5, MPFR_RNDD);
				mpfr_mul(t, t, s, MPFR_RNDU);	
			}	
     
			//mpfr_out_str (stdout, 10, 10, t, MPFR_RNDD);
			//putchar(' ');
			//mpfr_out_str(stdout, 10, 10, t_init, MPFR_RNDD);
			//mpfr_out_str (stdout, 10, 10, u, MPFR_RNDD);      
			//putchar (' ');
			mpfr_mul(t_final, t, t_init, MPFR_RNDD);
			//mpfr_out_str(stdout, 10, 10, t_final, MPFR_RNDD);
			//putchar(' ');
			m= mpfr_get_flt (t_final, MPFR_RNDD);
			//printf("%.10f\n",m);
			results[j] = results[j] + m;
		
			
				
			//mpfr_out_str (gnuplot, 10, 10, t, MPFR_RNDD);
			
		    
			//mpfr_set_d (u, 1.0, MPFR_RNDD);
	     
			//mpfr_add (s, s, u, MPFR_RNDD);
		}
		//mpfr_set_d(t_init, 0, MPFR_RNDD);
		//mpfr_set_d(t, 0, MPFR_RNDD);
	}
	printf("%.10f\n", xe[0]/numbi);
	for(int iter=0; iter<50; iter++){
		//temp = ((results[iter]/1000) - xs ) / ys;
		//temp=temp/10000;
		//temp = (temp-xs)/ys;
		//results[iter] = temp;
		fprintf(gnuplot_9u_cor, "%d %.10f\n", (iter+1), ((results[iter]/numbi)));// - ((xe[0]/numbi) * (xe[0]/numbi))));// / ((results[0]/numbi) - ((xe[0]/numbi) * (xe[0]/numbi))));
		fprintf(gnuplot_9u_cor, "\n");
	}
	//printf ("Sum is ");
	//mpfr_out_str (stdout, 10, 0, s, MPFR_RNDD);
	//putchar ('\n');
	//fprintf(gnuplot, "e\n");
	fflush(gnuplot_9u_cor);
	mpfr_clear (s);
	mpfr_clear (t);
	mpfr_clear (start);
	mpfr_clear (init);
	mpfr_clear (t_init);
	mpfr_clear (t_final);

	mpf_clear(x);
	mpf_clear(probm);

	//mpfr_clear (u);
	mpfr_free_cache ();
	//printf("0.99 ensemble done");
	return 0;
}