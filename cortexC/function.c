#include "protos.h"

/* --------------------- Evaluate --------------------- */

void RK4(float inj, float time, neuron * cell, int ind) {
	int i, j;
	float injTot;
	synapse * synAux;
	
	// Calculando k1
	for (i = 0; i < cell->numeq; i++)
		waux[i] = cell->w[i];

	injTot = inj - synCurrent(cell);
	(*cell->pfx)(injTot);
	for (i = 0; i < cell->numeq; i++)
		ks[i][0] = h * fbuf[i];
	
	synAux = cell->headSyn;
	for (i = 0; i < cell->numSyns; i++){
		for (j = 0; j < synAux->numeq; j++)
			waux[synAux->posEdo + j] = synAux->w[j];

		(*synAux->pfx)(synAux);
		for (j = 0; j < synAux->numeq; j++)
			ks[synAux->posEdo + j][0] = h * fbuf[synAux->posEdo + j];

		if (i < cell->numSyns - 1) 
			synAux = synAux->nextSyn;
	}
	
	// Calculando k2
	
	for (i = 0; i < cell->numeq; i++)
		waux[i] = cell->w[i] + ks[i][0]/2;
    
	injTot = inj - synCurrent(cell);
	(*cell->pfx)(injTot);
	for (i = 0; i < cell->numeq; i++)
		ks[i][1] = h * fbuf[i];
	
	synAux = cell->headSyn;
	for (i = 0; i < cell->numSyns; i++){
		for (j = 0; j < synAux->numeq; j++)
			waux[synAux->posEdo + j] = synAux->w[j] + ks[synAux->posEdo + j][0]/2;

		(*synAux->pfx)(synAux);
		for (j = 0; j < synAux->numeq; j++)
			ks[synAux->posEdo + j][1] = h * fbuf[synAux->posEdo + j];

		if (i < cell->numSyns - 1) 	
			synAux = synAux->nextSyn;
	}		
	
	// Calculando k3
	
	for (i = 0; i < cell->numeq; i++)
		waux[i] = cell->w[i] + ks[i][1]/2;

	injTot = inj - synCurrent(cell);
	(*cell->pfx)(injTot);
	for (i = 0; i < cell->numeq; i++)
		ks[i][2] = h * fbuf[i];
	
	synAux = cell->headSyn;
	for (i = 0; i < cell->numSyns; i++){
		for (j = 0; j < synAux->numeq; j++)
			waux[synAux->posEdo + j] = synAux->w[j] + ks[synAux->posEdo + j][1]/2;

		(*synAux->pfx)(synAux);
		for (j = 0; j < synAux->numeq; j++)
			ks[synAux->posEdo + j][2] = h * fbuf[synAux->posEdo + j];

		if (i < cell->numSyns - 1) 			
			synAux = synAux->nextSyn;
	}		
	
	// Calculando k4
	
	for (i = 0; i < cell->numeq; i++)
		waux[i] = cell->w[i] + ks[i][2];

	injTot = inj - synCurrent(cell);
	(*cell->pfx)(injTot);
	for (i = 0; i < cell->numeq; i++)
		ks[i][3] = h * fbuf[i];
	
	synAux = cell->headSyn;
	for (i = 0; i < cell->numSyns; i++){
		for (j = 0; j < synAux->numeq; j++)
			waux[synAux->posEdo + j] = synAux->w[j] + ks[synAux->posEdo + j][2];

		(*synAux->pfx)(synAux);
		for (j = 0; j < synAux->numeq; j++)
			ks[synAux->posEdo + j][3] = h * fbuf[synAux->posEdo + j];

		if (i < cell->numSyns - 1) 	
			synAux = synAux->nextSyn;
	}

	// Calculando w's novos
	
	for (i = 0; i < cell->numeq; i++)
		cell->w[i] += (ks[i][0] + 2*ks[i][1] + 2*ks[i][2] + ks[i][3])/6;

	(*cell->pcheck)(time, cell, ind);
	
	synAux = cell->headSyn;
	int posEdo;
	for (i = 0; i < cell->numSyns; i++){
		for (j = 0; j < synAux->numeq; j++){
			posEdo = synAux->posEdo + j;
			synAux->w[j] += (ks[posEdo][0] + 2*ks[posEdo][1] + 2*ks[posEdo][2] + ks[posEdo][3])/6;
		}
		checkEvtSyn(synAux,time);
		if (i < cell->numSyns - 1) 
			synAux = synAux->nextSyn;
	}	
}

/* --------------------- RsCell --------------------- */

void fxRsCell(float inj) {
	fbuf[0] = (0.7 * (waux[0] + 60.0) * (waux[0] + 40.0) - waux[1] + inj)/100.0;
	fbuf[1] = 0.03 * (-2.0 * (waux[1] + 60.0) - waux[1]);
}

void peakRsCell(float time, neuron * cell, int ind) {
    if (cell->w[0] >= 30.0) {
        cell->w[0] = -50.0;
        cell->w[1] = cell->w[1] + 100.0;
        recEvent(ind, time);
    }
} 

void createRsCell(neuron * cell){
	int i, j;
	cell->numeq = 2;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxRsCell;
	cell->pcheck = peakRsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -60.0;	cell->w[1] = 0;
}	

/* --------------------- IbsCell --------------------- */

void fxIbsCell(float inj) {
	fbuf[0] = (1.2 * (waux[0] + 75.0) * (waux[0] + 45.0) - waux[1] + inj)/150.0;
	fbuf[1] = 0.01 * (5.0 * (waux[0] + 75.0) - waux[1]);
}

void peakIbsCell(float time, neuron * cell, int ind) {
    if (cell->w[0] >= 30.0) {
        cell->w[0] = -56.0;
        cell->w[1] = cell->w[1] + 130.0;
        recEvent(ind, time);
    }
} 

void createIbsCell(neuron * cell){
	int i, j;
	cell->numeq = 2;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxIbsCell;
	cell->pcheck = peakIbsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -75.0;	cell->w[1] = 0;
}	
	
/* --------------------- ChsCell --------------------- */

void fxChsCell(float inj) {
	fbuf[0] = (1.5 * (waux[0] + 60.0) * (waux[0] + 40.0) - waux[1] + inj)/50.0;
	fbuf[1] = 0.03 * (1.0 * (waux[0] + 60.0) - waux[1]);
}

void peakChsCell(float time, neuron * cell, int ind) {
    if (cell->w[0] >= 30.0) {
        cell->w[0] = -40.0;
        cell->w[1] = cell->w[1] + 150.0;
        recEvent(ind, time);
    }
} 

void createChsCell(neuron * cell){
	int i, j;
	cell->numeq = 2;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxChsCell;
	cell->pcheck = peakChsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -60.0;	cell->w[1] = 0;
}

/* --------------------- FsCell --------------------- */

void fxFsCell(float inj) {
    fbuf[0] = ((waux[0] + 55.0) * (waux[0] + 40.0) - waux[1] + inj)/20.0;
    if (waux[0] < -55.0)
        fbuf[1] = - 0.2 * waux[1];
    else
        fbuf[1] = 0.2 * (0.025 * pow((waux[0] + 55.0),3) - waux[1]);
}

void peakFsCell(float time, neuron * cell, int ind) {
    if (cell->w[0] >= 25.0) {
        cell->w[0] = -45.0;
        recEvent(ind, time);
    }
} 

void createFsCell(neuron * cell){
	int i, j;
	cell->numeq = 2;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxFsCell;
	cell->pcheck = peakFsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -55.0;	cell->w[1] = 0;
}

/* --------------------- LtsCell --------------------- */

void fxLtsCell(float inj) {
    fbuf[0] = ((waux[0] + 56.0) * (waux[0] + 42.0) - waux[1] + inj)/100.0;
    fbuf[1] = 0.03 * (8.0 * (waux[0] + 56.0) - waux[1]);
}

void peakLtsCell(float time, neuron * cell, int ind) {
    if (cell->w[0] >= 40.0 - 0.1*cell->w[1]) {
        cell->w[0] = -53.0 + 0.04*cell->w[1];
        cell->w[1] = cell->w[1] + 20.0; 
        recEvent(ind, time);
    }
} 

void createLtsCell(neuron * cell){
	int i, j;
	cell->numeq = 2;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxLtsCell;
	cell->pcheck = peakLtsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -56.0;	cell->w[1] = 0;
}

/* --------------------- LsCell --------------------- */

void fxLsCell(float inj) {
	fbuf[0] = (0.3 * (waux[0] + 66.0) * (waux[0] + 40.0)  + 1.2*(waux[2] - waux[0]) - waux[1] + inj)/20.0;
    fbuf[1] = 0.17 * (5.0 * (waux[0] + 66.0) - waux[1]);
    fbuf[2] = (waux[0]-waux[2])/100.0;
}

void peakLsCell(float time, neuron * cell, int ind) {
	if (cell->w[0] >= 30.0) {
		cell->w[0] = -45.0;
		cell->w[1] = cell->w[1] + 100.0;
		recEvent(ind, time);
	}
} 

void createLsCell(neuron * cell){
	int i, j;
	cell->numeq = 3;
	cell->w = (double *) malloc(cell->numeq*sizeof (double));
	cell->pfx = fxLsCell;
	cell->pcheck = peakLsCell;
	cell->headSyn = NULL;
	cell->numSyns = 0;	
	cell->w[0] = -66.0;	cell->w[1] = 0; cell->w[2] = -66.0;
}

/* --------------------- Synapse --------------------- */

void createSyn(synapse * syn, short int type, short int delay, float weight, unsigned char shortTerm, unsigned char longTerm){
	int i, j;
	syn->numeq = 1;
	syn->weight = weight;
	syn->w = (double *) malloc(syn->numeq*sizeof (double));
	syn->delayType = (delay & 0xff) << 8 | (type & 0xff);
	syn->spkNum = 0;
	if (type == 0)
		syn->pfx = fxAmpa;
	if (type == 1)
		syn->pfx = fxNmda;
	if (type == 2)
		syn->pfx = fxGabaA;
	if (type == 3)
		syn->pfx = fxGabaB;
	syn->w[0] = 0.0;
	syn->preNeuron = -1;
	syn->nextSyn = NULL;
}

void fxGabaA(synapse * syn){
	int posEdo = syn->posEdo;
	fbuf[posEdo] = - (waux[posEdo]/6.0);
}

void fxGabaB(synapse * syn){
	int posEdo = syn->posEdo;
	fbuf[posEdo] = - (waux[posEdo]/150.0);
}

void fxAmpa(synapse * syn){
	int posEdo = syn->posEdo;
	fbuf[posEdo] = - (waux[posEdo]/5.0);
}

void fxNmda(synapse * syn){
	int posEdo = syn->posEdo;
	fbuf[posEdo] = - (waux[posEdo]/150); 
}

double synCurrent(neuron * cell) {
	double isyn = 0.0;
	int i = 0;
	int posEdo = cell->numeq;
	synapse * synAux = cell->headSyn;
	for (i = 0; i < cell->numSyns; i++){
		if ((synAux->delayType & 0xff) == 0)
			isyn += waux[posEdo]*(waux[0] - 0.0);	
		if ((synAux->delayType & 0xff) == 1)
			isyn += waux[posEdo]*pow((waux[0] + 80)/60,2)/(1 + pow((waux[0] + 80)/60,2))*(waux[0] - 0.0);
		if ((synAux->delayType & 0xff) == 2)
			isyn += waux[posEdo]*(waux[0] + 70.0);
		if ((synAux->delayType & 0xff) == 3)
			isyn += waux[posEdo]*(waux[0] + 90.0);	
		if (i < cell->numSyns - 1) 	
			synAux = synAux->nextSyn;
		posEdo += synAux->numeq;
	}		
	return isyn;
} 

void checkEvtSyn(synapse * syn, float time){
	float s;
	if (spkTime[syn->preNeuron][syn->spkNum] > 0.0){
		s = time - spkTime[syn->preNeuron][syn->spkNum] - ((syn->delayType & 0xff00) >> 8);
		if (s > 0){
			syn->w[0] += syn->weight;
			syn->spkNum++;
		}
	}        
	return;	
}

void recEvent(int ind, float time){
	float * tmp;
	int i;
	if (lengthSpk[ind] < maxSpkLength[ind]){
		spkTime[ind][lengthSpk[ind]] = time;
		lengthSpk[ind]++;
	}
	else{
		maxSpkLength[ind] += 50;
		tmp = (float *) malloc(maxSpkLength[ind] * sizeof( float ));
		for(i = 0; i < maxSpkLength[ind] - 50; i++){
			tmp[i] = spkTime[ind][i];
		}
		free(spkTime[ind]);
		spkTime[ind] = tmp;
		spkTime[ind][lengthSpk[ind]] = time;
	}
}

void createConnection(neuron * pos, synapse * syn, int pre){
	synapse * synTmp;
	int posicao = pos->numeq, j;
	if (pos->numSyns == 0){
		pos->headSyn = syn;
		syn->posEdo = posicao;
	}
	else{
		synTmp = pos->headSyn;
		posicao += synTmp->numeq;
		for (j = 1; j < pos->numSyns; j++){
			synTmp = synTmp->nextSyn;
			posicao += synTmp->numeq;
		}
		syn->posEdo = posicao;
		synTmp->nextSyn = syn;
	}
	pos->numSyns++;
	syn->preNeuron = pre;	
}

/* --------------------- Rede --------------------- */

void createRegNet (int nneuron, int numConn, int first, int nneuronNo, neuron * cells){    
    int conNeuron = (int) round(numConn/nneuron);
   	int k, i, j, pos, nEqMax = 10.000;
	synapse * syn;
	
	for (i = 0; i < nneuron; i++)
		createRsCell(&cells[i]);
	
	for (k = 0; k < nneuron; k++){
		for (i = 0; i < conNeuron; i++){
			pos = (i - conNeuron/2) + k;
			if (pos < 0)
				pos = nneuron + pos;
			if (pos > nneuron - 1)
				pos = pos - nneuron;
			
			syn = (synapse *) malloc(sizeof(synapse));
			createSyn(syn, 0, 1, 1.0, 0, 0);
			createConnection(&cells[pos], syn, k);
			
		}
	}
	
	ks = (double **) malloc( nEqMax * sizeof( double * ) );
	for(i = 0; i < nEqMax; i++)
		ks[i]= (double *) malloc( 4 * sizeof( double ) );		
	waux = (double *) malloc( nEqMax * sizeof( double ) );
	fbuf = (double *) malloc( nEqMax * sizeof( double ) );	
	
	return;	
}

void createCortexNet (int nneuron, float rescaleFac, float sizeNet, int * nLs, neuron * cells, int initNeuron, int endNeuron) {
    
    int i, j, k, l, p, count, aux;
    double scaleL1, scaleL23, scaleL4, scaleL5, scaleL6, scalePre, scalePos, sort;
    unsigned char type, typeCells[nneuron];
    synapse * syn;
    connPos * mtxPos = (connPos *) malloc(nneuron*sizeof(connPos));;
    
    //definindo escala por camada
    scaleL1 = sizeNet/(nLs[0] - 1);
    scaleL23 = sizeNet/(nLs[1] - 1);
    scaleL4 = sizeNet/(nLs[2] - 1);
    scaleL5 = sizeNet/(nLs[3] - 1);
    scaleL6 = sizeNet/(nLs[4] - 1);
    
    //carrega arquivos de conexao e arborizacao axonal
    FILE * connFile, * axonFile;
    double connDat[33][21], axonDat[17][6];    
    connFile = fopen("dataConn.dat","r");
    axonFile = fopen("dataAxon.dat","r");    
    for(i = 0; i < 33*21; i++){
    	fscanf(connFile,"%lf", &connDat[i/21][i%21]);
    }   
    for(i = 0; i < 17*6; i++){
    	fscanf(axonFile,"%lf", &axonDat[i/6][i%6]);
    	axonDat[i/6][i%6] *= rescaleFac;    	
    }       
    fclose(connFile);
    fclose(axonFile);
    
    //Definindo por camada, tipo de cada neurônio   
    count = 0;
    for (i = 0; i < nLs[0]*nLs[0]; i++) { 
    	typeCells[count] = 0;
	    mtxPos[count].length = 0;	
    	count++;
    }
    
    for (i = 0; i < nLs[1]*nLs[1]; i++) { 
    	sort = rand0();
    	if (sort < 0.78)
	    	typeCells[count] = 1;
    	else if (sort >= 0.78 && sort < 0.87)
	    	typeCells[count] = 2;    		
    	else
	    	typeCells[count] = 3;    		
    	mtxPos[count].length = 0;   	
    	count++;  	
    }
    
    for (i = 0; i < nLs[2]*nLs[2]; i++) { 
    	sort = rand0();
    	if (sort < 0.27)
	    	typeCells[count] = 4;
    	else if (sort >= 0.27 && sort < 0.54)
	    	typeCells[count] = 5;    		
    	else if (sort >= 0.54 && sort < 0.81)
	    	typeCells[count] = 6;    		
    	else if (sort >= 0.81 && sort < 0.96)
	    	typeCells[count] = 7;    		
    	else
	    	typeCells[count] = 8;    		
	    mtxPos[count].length = 0;	 	
    	count++;     	
    }
    
    for (i = 0; i < nLs[3]*nLs[3]; i++) {
    	sort = rand0();
    	if (sort < 0.64)
	    	typeCells[count] = 9;
    	else if (sort >= 0.64 && sort < 0.81)
	    	typeCells[count] = 10;    		
    	else if (sort >= 0.81 && sort < 0.89)
	    	typeCells[count] = 11;    		
    	else
	    	typeCells[count] = 12;    		
	    mtxPos[count].length = 0;	  	
    	count++;     	
    }
    
    for (i = 0; i < nLs[4]*nLs[4]; i++) {
    	sort = rand0();
    	if (sort < 0.62)
	    	typeCells[count] = 13;
    	else if (sort >= 0.62 && sort < 0.82)
	    	typeCells[count] = 14;	
    	else if (sort >= 0.82 && sort < 0.91)
	    	typeCells[count] = 15;	
    	else
	    	typeCells[count] = 16;
	    mtxPos[count].length = 0;
    	count++;     	
    }
    
	/* A info do arquivo de conexao é lida linha a linha (p). Para cada linha
	varre-se todos o neuronios, aqueles que pertencerem ao tipo especificado
	pela linha é analisadq a info da linha em questao. Com a leitura é coletado 
	o indice dos neuronios potencialmente pre-sinapticos (analisando o tipo e
	distância (arquivo de arborizacao axonal) ). Após isso é sorteado dessa lista
	a quantidade de conexoes especificadas no arquivo de conexao */    
    
	int xPre, yPre, xPos, yPos, num, typePre, typePos, beginI, endI;
	int beginK, endK, initPre, initPos, nlayerPre, nlayerPos, layerPre, layerPos;		
	double re, dist;
	int listNeuron[10000];
	int listDelay[10000];
	int numPre = 0;
	
	int * _preSyn;
	unsigned short int * _delayType;
	float * _weight;
	
	// varredura das linhas do arquivo de conexao
	for (p = 0; p < 33; p++){
		typePos = (int) connDat[p][0];
		layerPos = (int) connDat[p][1];	
		// varredura de todos o neuronios
		if (typePos == 0){
			beginI = 0;
			endI = nLs[0]*nLs[0];
		}
		else if (typePos == 1 || typePos == 2 || typePos == 3){
			beginI = nLs[0]*nLs[0];
			endI =  nLs[0]*nLs[0] + nLs[1]*nLs[1];				
		}
		else if (typePos == 4 || typePos == 5 || typePos == 6 || typePos == 7 || typePos == 8){
			beginI = nLs[0]*nLs[0] + nLs[1]*nLs[1];
			endI =  nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];		
		}
		else if (typePos == 9 || typePos == 10 || typePos == 11 || typePos == 12){
			beginI = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];	
			endI = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];		
		}
		else{
			beginI = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];	
			endI = nneuron;	
		}		
		for (i = beginI; i < endI; i++){
			//verifica se neuronio pos-sinaptico é o tipo especificado na linha
			if (typeCells[i] == typePos){
				// define a camada, escala e offset na lista de neuronios respectivo ao neuronio pos em questao			
				if (typePos == 0){
					initPos = 0;
					nlayerPos = nLs[0];
					scalePos = scaleL1;
				}
				else if (typePos == 1 || typePos == 2 || typePos == 3){
					initPos = nLs[0]*nLs[0];
					nlayerPos = nLs[1];	
					scalePos = scaleL23;				
				}
				else if (typePos == 4 || typePos == 5 || typePos == 6 || typePos == 7 || typePos == 8){
					initPos = nLs[0]*nLs[0] + nLs[1]*nLs[1];
					nlayerPos = nLs[2];	
					scalePos = scaleL4;				
				}
				else if (typePos == 9 || typePos == 10 || typePos == 11 || typePos == 12){
					initPos = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];
					nlayerPos = nLs[3];	
					scalePos = scaleL5;				
				}
				else{
					initPos = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];
					nlayerPos = nLs[4];		
					scalePos = scaleL6;		
				}			
				//define posicao (sem escala) do neuronio pos	
				xPos = (i - initPos)%nlayerPos;
				yPos = (i - initPos)/nlayerPos;	
				//varre lina do arquivo de conexao, coletando os neuronios pre-sinapticos de um determinado tipo no qual axonio atinge pos sinaptico
				for (j = 4; j < 21; j++){					
					if (connDat[p][j] > 0){
						//coleta informacao da arborizacao axonal do neuronio pre na camada informada pelo arquivo de conexao
						re = axonDat[j - 4][layerPos];
						//quanto neuronios pre de um detrminado tipo realizam sinapse neste neuronio pos
						count = (int) round(connDat[p][3]*(connDat[p][j]/(100.0*rescaleFac)));
						if (j == 4){
							beginK = 0;
							endK = nLs[0]*nLs[0];
						}
						else if (j == 5 || j == 6 || j == 7){
							beginK = nLs[0]*nLs[0];
							endK =  nLs[0]*nLs[0] + nLs[1]*nLs[1];				
						}
						else if (j == 8 || j == 9 || j == 10 || j == 11 || j == 12){
							beginK = nLs[0]*nLs[0] + nLs[1]*nLs[1];
							endK =  nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];		
						}
						else if (j == 13 || j == 14 || j == 15 || j == 16){
							beginK = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];	
							endK = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];		
						}
						else{
							beginK = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];	
							endK = nneuron;	
						}						
						for (k = beginK; k < endK; k++){	
							typePre = typeCells[k];	
							// define a camada, escala e offset na lista de neuronios respectivo ao neuronio pos em questao
							if(typePre == j - 4){
								if (typePre == 0){
									initPre = 0;
									nlayerPre = nLs[0];
									layerPre = 1;
									scalePre = scaleL1;
								}
								else if (typePre == 1 || typePre == 2 || typePre == 3){
									initPre = nLs[0]*nLs[0];
									nlayerPre = nLs[1];	
									layerPre = 2;
									scalePre = scaleL23;				
								}
								else if (typePre == 4 || typePre == 5 || typePre == 6 || typePre == 7 || typePre == 8){
									initPre = nLs[0]*nLs[0] + nLs[1]*nLs[1];
									nlayerPre = nLs[2];	
									layerPre = 3;	
									scalePre = scaleL4;			
								}
								else if (typePre == 9 || typePre == 10 || typePre == 11 || typePre == 12){
									initPre = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2];
									nlayerPre = nLs[3];		
									layerPre = 4;
									scalePre = scaleL5;			
								}
								else{
									initPre = nLs[0]*nLs[0] + nLs[1]*nLs[1] + nLs[2]*nLs[2] + nLs[3]*nLs[3];
									nlayerPre = nLs[4];		
									layerPre = 5;
									scalePre = scaleL6;		
								}	
								//define posicao (sem escala) do neuronio pre				
								xPre = (k - initPre)%nlayerPre;
								yPre = (k - initPre)/nlayerPre;	
								//verifica se neuronio pre atinge pos	
								dist = pow(xPre*scalePre-xPos*scalePos,2) + pow(yPre*scalePre-yPos*scalePos,2);
								if (dist < re*re){
									listNeuron[numPre] = k;
									dist = sqrt(dist + pow(abs(layerPre - layerPos)*0.25,2));	
									if (dist/0.1 < 1.0)
										listDelay[numPre] = 1; 
									else if (dist/0.1 > 20.0)
										listDelay[numPre] = 20;
									else
										listDelay[numPre] = (int) round(dist/0.1); 
									numPre++;
								}
							}
						}
						
						//coletados neuronios potencialmente pre-sinapticos é realizado sorteio 
						for (k = 0; k < count; k++){		
							sort = round(rand0()*(numPre-1));					
							num = listNeuron[(int) sort];
							
							aux = search(&mtxPos[num], i, listDelay[(int) sort]);
							
							if(aux < 0){
								mtxPos[i].length++;
								type = typeCells[i];
								_preSyn =  (int *) malloc(mtxPos[i].length*sizeof(int));
								_delayType = (unsigned short int *) malloc(mtxPos[i].length*sizeof(unsigned short int));
								_weight = (float *) malloc(mtxPos[i].length*sizeof(float));
								
								for (l = 0; l < (mtxPos[i].length - 1); l++){			
									_preSyn[l] = mtxPos[i].preSyn[l];
									_weight[l] = mtxPos[i].weight[l];
									_delayType[l] = mtxPos[i].delayType[l];									
								}
								
								_preSyn[l] = num;
								_weight[l] = 1.0;
								
								if (type == 0 || type == 2 || type == 3 || type == 7 || type == 8 || type == 11 || type == 12 || type == 15 || type == 16){
									if(rand0() < 0.5){
										_delayType[l] = (listDelay[(int) sort] & 0xff) << 8 | (3 & 0xff);
									}
									else{
										_delayType[l] = (listDelay[(int) sort] & 0xff) << 8 | (2 & 0xff);
									}								
								}
								else{
									if (rand0() < 0.27){
										_delayType[l] = (listDelay[(int) sort] & 0xff) << 8 | (1 & 0xff);
									}
									else{
										_delayType[l] = (listDelay[(int) sort] & 0xff) << 8 | (0 & 0xff);
									}
								}
								
								free(mtxPos[i].preSyn);
								free(mtxPos[i].weight);
								free(mtxPos[i].delayType);
								mtxPos[i].preSyn = _preSyn;
								mtxPos[i].weight = _weight;
								mtxPos[i].delayType = _delayType;
							}
							else
								mtxPos[i].weight[aux]++;				
						}
						numPre = 0;
					}					
				}	
			}
		}	
	}

	for (i = initNeuron; i <= endNeuron; i++){
		type = typeCells[i];
		if (type == 1 || type == 4 || type == 5 || type == 6 || type == 9 || type == 10 || type == 13 || type == 14)
			createRsCell(&cells[i - initNeuron]);
		else if (type == 2 || type == 7 || type == 11 || type == 15)
			createFsCell(&cells[i - initNeuron]);
		else if (type == 3 || type == 8 || type == 12 || type == 16)
			createLtsCell(&cells[i - initNeuron]);
		else
			createLsCell(&cells[i - initNeuron]);
		
		for(j = 0; j < mtxPos[i].length; j++){
			syn = (synapse *) malloc(sizeof(synapse));
			createSyn(syn, (mtxPos[i].delayType[j] & 0xff), (mtxPos[i].delayType[j] & 0xff00) >> 8, mtxPos[i].weight[j], 0, 0);
			createConnection(&cells[i-initNeuron], syn, mtxPos[i].preSyn[j]);			
		}
	}

	for(i = 0; i < nneuron; i++){
		free(mtxPos[i].preSyn);
		free(mtxPos[i].weight);
		free(mtxPos[i].delayType);	
	}
	free(mtxPos);
	
	return;
} 

int search(connPos * mtxPos, int num, int delay){
	for (int i = 0; i < mtxPos->length; i++){
		if(mtxPos->preSyn[i] == num && ((mtxPos->delayType[i]  & 0xff) << 8) == delay)
			return i;
	}
	return -1;
}

void setSeed(long value){
	idum = value;
}
    
float rand0()
{
	int j; 
	long k; 
	static long idum2=123456789; 
	static long iy=0; 
	static long iv[NTAB]; 
	float temp;
	if (idum <= 0) { 
		if (-(idum) < 1) idum=1; 
		else idum = -(idum); 
		idum2=(idum); 
		for (j=NTAB+7;j>=0;j--) {
			k=(idum)/IQ1; 
			idum=IA1*(idum-k*IQ1)-k*IR1; 
			if (idum < 0) idum += IM1; 
			if (j < NTAB) iv[j] = idum;
		} 
		iy=iv[0];
	} 
	k=(idum)/IQ1; 
	idum=IA1*(idum-k*IQ1)-k*IR1; 
	if (idum < 0) idum += IM1; 
	k=idum2/IQ2; 
	idum2=IA2*(idum2-k*IQ2)-k*IR2; 
	if (idum2 < 0) idum2 += IM2; 
	j=iy/NDIV; iy=iv[j]-idum2; iv[j] = idum; 
	if (iy < 1) iy += IMM1; 
	if ((temp=AM*iy) > RNMX) return RNMX; 
	else return temp;
}