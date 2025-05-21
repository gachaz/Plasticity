ACCESS INFORMATION
1. Licenses/restrictions placed on the data or code
	GNU Lesser General Public License

SUMMARY
This code was written in 2024 to generate all simulations from the article by Lambert, Achaz, Le Rouzic and Loison (2025)
"The Baldwin effect reloaded: Intermediate levels of phenotypic plasticity favor evolutionary rescue"
It is an individual-based birth-death simulator in continuous time. Using a standard Gillespie algorithm
to generate waiting times until the next event (that is birth or death). Mutations are only produced at birth.
individuals with  'b' genotypes can develop into BETA phenotypes, whereas 'a' can develop into 'BETA' in the second
environment.

All comments on this code should be addressed to Guillaume Achaz <guillaume.achaz@college-de-france.fr>


DATA & CODE FILE OVERVIEW
This data repository consist of one C source files and this README document, with the following code file
Code scripts and workflow
    1. plasticity_be.c
		to compile the C code, use cc (or gcc), typically
			gcc -o plasticity_be -O2 -Wall plasticity_be.c
		to generate a binary entitled 'plasticity_be'
		
		General usage is './plasticity_be [options]'

		To get the description of the options, type './plasticity_be -h' (also listed below)
		
		[general options]  
			-h   : print help and exit
			-x   : set random seed (otherwise use time()) for pseudo-random number generator
			-r # : number of replicates (def: 1)
			-S # : sample once at generation # after reaching K* (def: none)
			-s # : sample every # generations after reaching K* (def: none)  - also used to monitor the population

		[models options]  
			-K # : set carying capacity (def. 1000)
			-T # : set time shift (def. -1 = none). When T>0, it sets when the env changes from 0 to 1.
			-n # : set the initial # of a individuals at t=0 (def. K*, where K* is the eq value in env 0, see article)
			-N # : set the initial # of b individuals at t=0 (def. 0)
			-p # : plasticity to be a->BETA (def. 0.1). This the probability that an "A" genotype develops into a "BETA" phenotype in env 1.
			-g # : set gmax (max # of evts). This can be used to limit the simulation time.
			-G # : set tmax (max time). This can be used to limit the simulation time.

		[rate options]  
			-d # : set death rate of type A (def. 1.0)
			-D # : --- ----- ---- -- ---- B (def. 1.0)
			-b # : set birth rate of type A (def. 2.0 in env A, 0.5 in env B)
			-B # : --- ----- ---- -- ---- B (def. 0.5 -- --- A, 2.0 -- --- B)
			-c # : set competition rate of type A (def. 1.0)
			-C # : --- ----------- ---- -- ---- B (def. 1.0)
			-m # : set mutation rate of a2b (mut rate = #/K, def. 1.0 that is 1 mut/gen)
			-M # : --- -------- ---- -- b2a (def. 0.0)


SOFTWARE VERSIONS
	C standard stand-alone (no specific library is required)