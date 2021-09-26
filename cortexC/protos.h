#ifndef _PROTOS_
#define _PROTOS_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define IM1 2147483563 
#define IM2 2147483399 
#define AM (1.0/IM1) 
#define IMM1 (IM1-1) 
#define IA1 40014 
#define IA2 40692 
#define IQ1 53668 
#define IQ2 52774 
#define IR1 12211 
#define IR2 3791 
#define NTAB 32
#define NDIV (1+IMM1/NTAB) 
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

typedef struct synapse{
	double * w;
	unsigned short int numeq;
	unsigned short int delayType;
	unsigned short int posEdo;
	unsigned short int spkNum;
	float weight;
	struct synapse * nextSyn;
	int preNeuron;
	void (*pfx) (struct synapse * syn);
} synapse;

typedef struct neuron{
	double * w;
	short int numeq;
	synapse * headSyn;
	short int numSyns;
	void (*pfx) (float inj);
	void (*pcheck) (float time, struct neuron * cell, int ind);
} neuron;

typedef struct connPos{
	int * preSyn;
	unsigned short int * delayType;
	float * weight;	
	unsigned short int length;
} connPos;

void RK4(float inj, float time, neuron * cell, int ind);
void createRsCell(neuron * cell);
void fxRsCell(float inj);
void peakRsCell(float time, neuron * cell, int ind);
void createChsCell(neuron * cell);
void fxChsCell(float inj);
void peakChsCell(float time, neuron * cell, int ind);
void createIbsCell(neuron * cell);
void fxIbsCell(float inj);
void peakIbsCell(float time, neuron * cell, int ind);
void createFsCell(neuron * cell);
void fxFsCell(float inj);
void peakFsCell(float time, neuron * cell, int ind);
void createLtsCell(neuron * cell);
void fxLtsCell(float inj);
void peakLtsCell(float time, neuron * cell, int ind);
void createLsCell(neuron * cell);
void fxLsCell(float inj);
void peakLsCell(float time, neuron * cell, int ind);

double synCurrent(neuron * cell);
void recEvent(int ind, float time);
void checkEvtSyn(synapse * syn, float time);

void createSyn(synapse * syn, short int type, short int delay, float weight, unsigned char shortTerm, unsigned char longTerm);
void fxGabaA(synapse * syn);
void fxGabaB(synapse * syn);
void fxAmpa(synapse * syn);
void fxNmda(synapse * syn);
void createConnection(neuron * pos, synapse * syn, int pre);

void createRegNet (int nneuron, int numConn, int first, int nneuronNo, neuron * cells);
void createCortexNet (int nneuron, float rescaleFac, float sizeNet, int * nLs, neuron * cells, int initNeuron, int endNeuron);
int search(connPos * mtxPos, int num, int delay);

void setSeed(long value);
float rand0();

extern float ** spkTime;
extern float ** spkBuffer;
extern unsigned short int * lengthSpk;
extern unsigned short int * maxSpkLength;
extern double ** ks;
extern double * waux;
extern double * fbuf;
extern double h;
extern long int idum;

#endif