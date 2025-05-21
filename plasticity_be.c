/*
        Copyright (C) 2024 to Guillaume Achaz.

        This program is free software; you can redistribute it and/or
        modify it under the terms of the GNU Lesser General Public License
        as published by the Free Software Foundation; either version 2.1
        of the License, or (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU Lesser General Public License for more details.

        You should have received a copy of the GNU Lesser General Public License
        along with this program; if not, write to the Free Software
        Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
        for more information, please contact Guillaume Achaz <guillaume.achaz@college-de-france.fr>
*/


/*
	This code was written by G Achaz in 2022-2024
	It performs individual birth-death population simulations
	to test for the effect of plasticity.
	
	It has be commented and streamlined so that usage and code reading
	should be easy. As a stand-alons, it has no dependance on any library.
	
	The corresponding ms is
	"The Baldwin effect reloaded: Intermediate levels of
      phenotypic plasticity favor evolutionary rescue"
      
    authors: Lambert, Achaz, Le Rouzic and Loison
                              
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <unistd.h>

#define HISTOSIZE 1000     // no more than HISTOSIZE 'b' genotypes inv in 0 -- for mut-sel balance.
#define MAXTYPE 10000      // no more than MAXTYPE of geno b are allowed

char verbose;

/*
	A structure to follow segregating types in the 'b' pop (and their timing of appearance)
*/
struct all_types {
	int histo[MAXTYPE];     // how many of each descendants of each type among the 'b' genotypes
	int min;                // there are currently no type below this one
	int max;                //   "    "      "     "    "  above   "   "
	float origin[MAXTYPE];  // timing of first appearance of each mutant
};


/*
	Output all the geno/PHENO
*/
void PrintPop( int g, float time, int *pop ){
	printf("g: %d t: %f pop: %d %d %d %d\n", g, time, pop[0], pop[1], pop[2], pop[3]);
}
/*
	Output the b types segregating
*/
void PrintTypes( struct all_types *my_types, int geno_b ){
	
	int i, sum_b=0;
	
	printf("[min, max]: [ %d , %d ]\n", my_types->min , my_types->max );
	
	for(i=my_types->min; i<=my_types->max; i++)
	{
		printf("histo[%d]: %d ; t_mut: %f\n", i, my_types->histo[i], my_types->origin[i] );
		sum_b+=my_types->histo[i];
	}
	printf("sum: %d ; genob: %d\n", sum_b, geno_b);
}

/***
	Random number generators
***/

/*
	uni_float in [0,1[
	!! checked on MACOSX !!
*/
double unirand( void ){

	double r=(double)random();
	return r/pow(2.0,31);
}

/*
	Return a random uniform in [min,max]
*/
int unirand_int( int min, int max){
	return min + (int) floor( unirand()*(max-min+1.0) );
}

/*
	exp_float in [0,\infty)
*/
double exprand( float rate ){
	return -log( 1-unirand( ) )/rate ;
}



/*
	This draws a random type among all 'b' individuals with
	respect to their abundances
*/
int pick_type_of_b( struct all_types *my_types, int number_of_b ){

	int t = my_types->min,                    // dummy counter iterating over all types of genob
	    r = unirand_int(0,number_of_b-1);     // random int in [0,number_of_b[
	
	while ( r >= my_types->histo[t] && t <= my_types->max ){
		r-=my_types->histo[t];
		t++;
	}

	return t;
}

float get_rescue_time( struct all_types my_types ){

	int t,
	    max_type = my_types.min;
	
	for (t = my_types.min+1 ; t <= my_types.max ; t++)
		if( my_types.histo[t] > my_types.histo[max_type])
			max_type = t;

	return my_types.origin[max_type];
}



/*
	pop is {geno_a>pheno_ALPHA, geno_a>pheno_BETA, b2ALPHA, b2BETA}
	float P_a2BETA is the current probability of having an a2BETA (0 before env. shift)
	Add an individual
*/
void Birth( int *pop, struct all_types *my_types, float *b, float P_a2BETA, float Q_b2ALPHA, float *m, float time ){
	
	float rate = (pop[0]+pop[2])*b[0]+(pop[1]+pop[3])*b[1];
	
	double r1=unirand(),   // select PHENO
	       r2=unirand(),   // check for mutations
	       r3=unirand();   // once geno is chosen, select PHENO
	
	int t;


	if ( (r1 < (pop[0]*b[0]+pop[1]*b[1])/rate && r2 >= m[0]) ||        // parent is 'a' + no mut.
	     (r1 >= (pop[0]*b[0]+pop[1]*b[1])/rate && r2 < m[1] )   )      // parent is 'b' + mut.
	{
		pop[ 0+( (r3<P_a2BETA)?1:0 ) ]++;                           // in both cases; the newborn is geno_a)
	}
	else
	{
		pop[ 2+( (r3<Q_b2ALPHA)?0:1 ) ]++;                           // otherwise, the newborn is geno_b

		/*
			Then, update Types
		*/
		if (r1 < (pop[0]*b[0]+pop[1]*b[1])/rate && r2 < m[0] )   // parent is 'a' + mut,
		{
		
			my_types->max++;
			my_types->histo[my_types->max]=1;
			my_types->origin[my_types->max]=time;             // this is to follow when geno b first appear
			
		}
		else                                                    // parent is 'b' + no mut.
		{
		
			t = pick_type_of_b( my_types, pop[2]+pop[3]-1 ); // as the birth event is already in pop
			my_types->histo[ t ]++;
		}
		
	}
	
	return;
}


/*
	pop is {geno_a>pheno_ALPHA, geno_a>pheno_BETA, b2ALPHA, b2BETA}
	Kill an individual, if genob, update Types
*/
void Death( int *pop, struct all_types *my_types, float *death_r ){
	
	float rate = (pop[0]+pop[2])*death_r[0]+(pop[1]+pop[3])*death_r[1];
	
	double r1=unirand(),   // select PHENO
	       r2=unirand();   // select geno
	
	int type,              // type of genob
	    flag_geno_b=0;     // becomes 1 when a geno b is killed
	

	if ( r1 < (pop[0]+pop[2])*death_r[0]/rate )     // kill an ALPHA phenotype
	{
	
		if( r2< pop[0]/(0.0+pop[0]+pop[2]) )   // kill an a2ALPHA
			pop[0]--;    
		else                                   // kill a b2ALPHA
		{
			pop[2]--;     
			flag_geno_b=1;
		}
	
	}
	else  // kill a BETA phenotype
	{
		if(  r2<pop[1]/(0.0+pop[1]+pop[3]) )  // kill an a2BETA
			pop[ 1 ]--;
		else
		{
			pop[ 3 ]--;                    // kill a b2BETA
			flag_geno_b=1;
		}
	
	}
	
	if(flag_geno_b)    // if a 'b' genotype is killed, update the types
	{
		type = pick_type_of_b( my_types, pop[2]+pop[3] );              // update types (a random one dies)
		my_types->histo[ type ]--;

		while( my_types->histo[my_types->min] == 0 && my_types->min<my_types->max)   // eventually update min
			my_types->min++;
		
		while( my_types->histo[my_types->max] == 0 && my_types->max>=my_types->min )   // eventually update max
			my_types->max--;

	}
	
	return;

}



/*
	Pick randomnly and implement a birth or a death, according to the rates
	return time to first event. Gillespie-like algorithm.
	pop is [ a2ALPHA, a2BETA, b2ALPHA, b2BETA ]
	my_types is used to follow the different b mutants
	b : birth rates of ALPHA and BETA ; d death rates ; c competition terms ; mutation prob
	p: proba that a genotype a->BETA
	q: proba that a genotype b->ALPHA (0 in the article)
*/
float OneTimeStep( int *pop, struct all_types *my_types, float *b, float *d, float *c,  float *m, float K, float p, float q, float time  ){

	int A = pop[0]+pop[2],   // total number of ALPHA
	    B = pop[1]+pop[3];   //  "      "       BETA
	    
	float death_r[2] = { d[0]+c[0]*(A+B-1.0)/K , d[1]+c[1]*(A+B-1.0)/K };
	
	float total_rate = A*b[0]+B*b[1]+A*death_r[0]+B*death_r[1] ;

	double r=unirand();
	
	if(verbose == 2)
		printf("r= %f ; birth= %f ; death= %f %f\n", r, (A*b[0]+B*b[1]), total_rate-(A*b[0]+B*b[1]), A*death_r[0]+B*death_r[1]);

	if (r < (A*death_r[0]+B*death_r[1])/total_rate)
	{
		//printf("death\n");
		Death( pop, my_types, death_r );
	}
	else
	{
		//printf("birth\n");
		Birth( pop, my_types, b, p, q, m, time );
	}
	
	
	return exprand(total_rate);
}



void PrintHelp( char *prog ){

	putchar('\n');
	printf("  Usage is '%s [options]'\n", prog);

	putchar('\n');
	printf("  [general options]  \n");
	printf("\t-h   : print this help and exit\n");
	printf("\t-x   : set random seed (otherwise use time())\n");
	printf("\t-r # : number of replicates (def: 1)\n");
	printf("\t-S # : sample once at generation # after reaching K* (def: none)\n");
	printf("\t-s # : sample every # generations after reaching K* (def: none)\n");

	putchar('\n');
	printf("  [models options]  \n");
	printf("\t-K # : set carying capacity (def. 1000)\n");
	printf("\t-T # : set time shift (def. -1 = none)\n");
	printf("\t-n # : set the initial # of a individuals at t=0 (def. K*)\n");
	printf("\t-N # : set the initial # of b individuals at t=0 (def. 0)\n");
	printf("\t-p # : plasticity to be a->BETA (def. 0.1)\n");
	printf("\t-g # : set gmax (max # of evts)\n");
	printf("\t-G # : set tmax (max time)\n");

	putchar('\n');
	printf("  [rate options]  \n");
	printf("\t-d # : set death rate of type A (def. 1.0)\n");
	printf("\t-D # : --- ----- ---- -- ---- B (def. 1.0)\n");
	printf("\t-b # : set birth rate of type A (def. 2.0 in env A, 0.5 in env B)\n");
	printf("\t-B # : --- ----- ---- -- ---- B (def. 0.5 -- --- A, 2.0 -- --- B)\n");
	printf("\t-c # : set competition rate of type A (def. 1.0)\n");
	printf("\t-C # : --- ----------- ---- -- ---- B (def. 1.0)\n");
	printf("\t-m # : set mutation rate of a2b (mut rate = #/K, def. 1.0)\n");
	printf("\t-M # : --- -------- ---- -- b2a (def. 0.0)\n");


}


int main(int argc, char ** argv)
{
	extern char verbose;
	extern char *optarg;
	extern int optind; 

	float K=1000;                   // sort of carrying capacity, but not exactly, as N -> K*(b-d)/c
	int gmax=-1;                    //
	float tmax=-1;                  // maximal #events OR time for simulations

	int rep=1;

	float T_shift=-1;               // time after which environment is shifting from 0 to 1
	char shift=0;                   // state variable, did it shift ?

	char report_time_rescue=0;      // output time to rescue.
	
	char sample_once=0,
	     perform_sampling=0;
	int sample_b=0;


	float b[2] = {2.0, 0.5};   // birth rates of {ALPHA, BETA} phenotypes, regardless of their genotype
	float d[2] = {1.0, 1.0};   // death rates
	float c[2] = {1.0, 1.0};   // competition terms rates


	float mu=1;
	float m[2] = {mu, 0};      // mutation prob: mut_a2b, mut_b2a  (this concerns genotypes !!)
	                           // !! it will be divided by K below !!

	int histo_b[HISTOSIZE]={0,};   // monitor the number of 'b' in the least favorable environnement.

	float P_a2BETA = 0.1;     // plasticity prob that genotype 'a' develops into a 'BETA' phenotype
	                          // in the second env. This is the 'adaptive plasticity'.
	
	int initial_a=0;
	int initial_b=0;
	int g0;
	      
	int seed=time(NULL);
	
	struct all_types my_types;  // this record all types currently in pop and their time of appearance

	verbose=0;

	int o;
	
	while ((o = getopt(argc, argv, "vVr:g:b:B:d:D:c:C:m:M:n:N:T:p:s:S:x:hK:tG:")) != -1) {
	
		switch( o ){
		
			case 'h':
				PrintHelp(argv[0]);
				exit(0);
		
			case 'v':
				verbose = 1;
				break;
			case 'V':
				verbose = 2;
				break;

			case 'r':
				rep = atoi(optarg);
				break;

			case 'g':
				gmax = atoi(optarg);
				break;

			case 'G':
				tmax = atof(optarg);
				break;

			case 'p':
				P_a2BETA = atof(optarg);
				break;


			case 'K':
				K = atof(optarg);
				break;


			/*
				Rates
			*/
			case 'b':
				b[0] = atof(optarg);
				break;

			case 'B':
				b[1] = atof(optarg);
				break;

			case 'c':
				c[0] = atof(optarg);
				break;

			case 'C':
				c[1] = atof(optarg);
				break;

			case 'd':
				d[0] = atof(optarg);
				break;

			case 'D':
				d[1] = atof(optarg);
				break;

			case 'm':
				m[0] = atof(optarg);
				break;

			case 'M':
				m[1] = atof(optarg);
				break;


			/*
				Environment
			*/
			case 'T':
				T_shift = atof(optarg);
				break;

			/*
				Sampling
			*/
			case 'S':
				sample_once = 1;
			case 's':
				sample_b   = atoi(optarg);
				perform_sampling = 1;
				break;

			case 'x':
				seed = atoi(optarg);
				break;

			case 'n':
				initial_a = atoi(optarg);
				break;

			case 'N':
				initial_b = atoi(optarg);
				break;

			case 't':
				report_time_rescue =1;
				break;

			default:
				printf("unknown option: %c, bye\n",(char)o),exit(1);
		
		}
	}
	
	argc -= optind;
	argv += optind;
	
	m[0]/=K;          // rescale the mutation rate by K
	m[1]/=K;
	
	/*
		Header
	*/
	printf("/*\n");
	printf("   [ sims ]\n");
	printf("   seed = %d\n", seed);
	if(gmax>=0)
		printf("   gmax = %d\n", gmax);
	if(tmax>=0)
		printf("   tmax = %.1f\n", tmax);
		
	printf("   rep = %d\n", rep);
	printf("   verbose = %d\n", verbose);
	if(sample_b>0){
		printf("   sample every: %d\n", sample_b);
		if(sample_once == 0)
			printf("   => sampling size: 10^%.1f\n", log(gmax/sample_b)/(log(10)));
		else
			printf("   => sampling size: 10^%.1f\n", log(rep)/(log(10)));
	}
	
	/*
		If there no 'a' genotypes at the begining, set them to the carrying capacity
	*/
	if(initial_a == 0)
		initial_a=(K*(b[0]-d[0])/c[0]);

	/*
		Header, second part
	*/
	printf("   [ parameters ]\n");
	printf("                  [A]   [B]\n");
	printf("   birth rates : %.3f %.3f\n", b[0], b[1]);
	printf("   death rates : %.3f %.3f\n", d[0], d[1]);
	printf("   compet terms: %.3f %.3f\n", c[0], c[1]);
	printf("   mutat. probs: %.3f %.3f\n", K*m[0], K*m[1]);
	printf("       [pop]\n");
	
	float xeq=(b[0]-d[0])/c[0];
	
	printf("   K  /  K*          : %.1f / %.1f\n", K, K*xeq);
	printf("   #geno_a_t0  : %d\n", initial_a);
	printf("   #geno_b_t0  : %d\n", initial_b);
	printf("   P(a->BETA)  : %f (only in environment 2)\n", P_a2BETA);
	
	if(T_shift>=0){
	printf("   T_shift     : %f\n", T_shift);
	printf("   (after T_shift birth rates of ALPHA and BETA are swapped)\n");
	}
	
	printf("*/\n");
	
	
	srandom( seed );

	int ext=0;
	float p=0;     // the plasticity, in practice: 0 in env-1 and P_a2BETA in env-2
	int fix_b=0;
	int flag=0;
	float t_resc=-1;

	for( int r=0; r<rep; r++ )
	{
	
	
		int g=0;
		float t=0;
		
		flag=0;
		t_resc=-1;
		            
		int pop[4] = { (int)((1.0-p)*initial_a), (int)(p*initial_a), 0, initial_b };  // { a2ALPHA, a2BETA, b2ALPHA, b2BETA }

		for(int i=0;i<MAXTYPE;i++) my_types.histo[i] = my_types.origin[i] = 0;  // reset the type of geno_b
		my_types.min=my_types.max=0;
		
		if( initial_b > 0)
			my_types.histo[0] = initial_b;   // {all initial b's carry type 0}
		else
			my_types.max=-1;
		
		if(verbose)PrintPop( g, t, pop );
		
		
		if(r%100 ==0)
			fprintf(stderr, "%d/%d\r", r, rep);
					
		
		/*
			While pop is not extinct
		*/
		while( pop[0]+pop[1]+pop[2]+pop[3]>0 ){

		
			/*
				Whenever, there is a gmax/tmax, check it
			*/
			if(gmax>0 && g>gmax)
				break;

			if(tmax>0 && t>tmax)
				break;


			/*
				Once env. is shifted after T_shift>=0,
				swap the birth rates of both types
			*/
			if(T_shift>=0 && t>T_shift)
			{	
				float tmp=b[0];
				b[0]=b[1];
				b[1]=tmp;
				T_shift=-1;
				shift=1;
				p=P_a2BETA;     // now a genotypes are plastic
			}
			
			/*
				Stop when the b genotypes have invaded
				after the environnemental shift
			*/
			if( shift && pop[2]+pop[3] >= K*(b[1]-d[1])/c[1] && flag==0)
			{
				fix_b++;    // count runs where it occured
				flag=1;
				t_resc = get_rescue_time( my_types );
				//break;
			}
			
			
			/*
				Usefull for sampling the b before the change of env
				to measure mut-sel balance
				a. perform_sampling = 1 : samping is ongoing, phase 1
				b. shift = 0 : the env did not change yet
				c. pop_0+pop_1>=equilibrium value. it is at steady state.
			*/
			if( perform_sampling==1 && shift==0 && pop[0]+pop[1]>=(int)(K*(b[0]-d[0])/c[0]) )
			{
				perform_sampling=2;
				g0=g;
			}
			
			/*
				a. Sampling is now ongoing (perform_sampling in phase 2)
				b. (g-g0)%sample_b sample every sample_b from g0
				c. shift=0, still in env 0.
			*/
			if(perform_sampling==2 && g>=sample_b && (g-g0)%sample_b==0 && shift == 0 ){
			
				if(pop[2]+pop[3]>=HISTOSIZE)
					fprintf(stderr, "cannot store #b in hist, bye\n"), exit(1);
//				printf("sample at gen %d\n",g);
			
				histo_b[ pop[2]+pop[3] ]++;
				
				if(sample_once == 1)       // in this sampling scheme, b are sample only once for each replicate.
				{
					perform_sampling=1;
					break;
				}
			}


			/*
				Time flies
			*/
			t += OneTimeStep( pop, &my_types, b, d, c, m, K, p, 0, t  );
			
			
			/*
				Sanity check
			*/
			if(pop[0]<0 || pop[1]<0 || pop[2]<0 || pop[3]<0)
			{
				printf("error at rep %d, gen:  %d\n",r, g);
				PrintPop( g, t, pop );
				exit(5);
			}

			
			
			if(verbose)
				PrintPop( g, t, pop );


			g++;


		}

		if( pop[0]+pop[1]+pop[2]+pop[3]==0 )
			ext++;
			
		if(report_time_rescue)
		{
			if(tmax<=0){
				printf("please select a tmax using -G option, bye\n");
				exit(4);
			}
		
			if( pop[0]+pop[1]+pop[2]+pop[3]==0 )
				printf("t_resc= %f (ext)\n", tmax*1.1);
			else
			{
				if( t>=tmax && flag==0 )
					printf("t_resc= %f (tmax)\n", tmax);
				else
					printf("t_resc= %f\n",t_resc);
			}
		}
		
//		printf("%d / %d \n",r,rep);

//		PrintPop( g, t, pop );
	}
	
	
	if(sample_once || sample_b>0){
		int i;
		int s=0;
		float theta = m[0]*b[0]*xeq*K/b[1],
		       rho  = b[1]/(d[1] + c[1]*xeq),
		       P;

		for(i=0;i<HISTOSIZE;i++)
			s+= histo_b[i];
			
		printf("Sample size=%d\n", s);
			
		for(i=0;i<20;i++)
		{
			P = pow(1-rho, theta) * pow(rho,i) * tgamma(i+theta)/(tgamma(i+1.0)*tgamma(theta));
			printf("%d %f %f\n",i,histo_b[i]/(s+0.0), P);	
		}
	}
	else
		printf("P_ext= %f [ %d ]; tmax: %f [ %d ];  bfix: %f [ %d ])\n",
		         ext/(double)rep, ext, (rep-ext-fix_b)/(double)rep, rep-ext-fix_b, fix_b/(double)rep, fix_b );

	return 0;
}