#include "protos.h"

float ** spkTime;
float ** spkBuffer;
unsigned short int * lengthSpk;
unsigned short int * maxSpkLength;
double ** ks;
double * waux;
double * fbuf;
double h = 0.01;
long int idum = -1089432789;

int main(int argc, char * argv[]) {

	float t = 0.0f, tmax = 1.0f;
	float rescaleFac = 100.0, sizeNet = 2.0;
	int Nstep = (int) (tmax - t)/h + 1, nNeuron = 3129;
	int i, j, nLs[5], maxSize = 1000;
	
	setSeed(-1089432789);
	
    nLs[0] = (int) round(sqrt(0.015*nNeuron));
    nLs[1] = (int) round(sqrt(0.337*nNeuron));
    nLs[2] = (int) round(sqrt(0.349*nNeuron));
    nLs[3] = (int) round(sqrt(0.076*nNeuron));
    nLs[4] = (int) round(sqrt(0.223*nNeuron));		
    
    nNeuron = 0;
    for (i = 0; i < 5; i++)
    	nNeuron += nLs[i]*nLs[i];

	neuron * cells = (neuron *) malloc(nNeuron * sizeof( neuron ));
	createCortexNet (nNeuron, rescaleFac, sizeNet, nLs, cells, 0, nNeuron-1);			

	spkTime = (float **) malloc(nNeuron * sizeof( float * ));
	for (i = 0; i < nNeuron; i++){
		spkTime[i] = (float *) malloc(50 * sizeof( float ));
	}
	
	for(j = 0; j < nNeuron; j++){
		for(i = 0; i < 50; i++)
			spkTime[j][i] = -1.0f;
	}
	
	ks = (double **) malloc( maxSize * sizeof( double * ) );
	for(i = 0; i < maxSize; i++)
		ks[i]= (double *) malloc( 4 * sizeof( double ) );		
	waux = (double *) malloc( maxSize * sizeof( double ) );
	fbuf = (double *) malloc( maxSize * sizeof( double ) );
	
	lengthSpk = (unsigned short int *) malloc(nNeuron * sizeof(unsigned short int));
	maxSpkLength = (unsigned short int *) malloc(nNeuron * sizeof(unsigned short int));
	memset(lengthSpk, 0, nNeuron*sizeof(unsigned short int));
	for(i = 0; i < nNeuron; i++)
		maxSpkLength[i] = 50;
		
	FILE * pfile = fopen("saida.dat", "w");
	for(i = 0; i < Nstep; i++){
		t += h;
		for(j = 0; j < nNeuron; j++)
			RK4(250.0, t, &cells[j], j);
		fprintf(pfile, "%f\t%lf\t%lf\n", t, cells[0].w[0], cells[1].w[0]);		
	}
	fclose(pfile);	
		
	return 0;
}