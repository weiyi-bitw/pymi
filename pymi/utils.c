/*
 * =====================================================================================
 *
 *       Filename:  utils.c
 *
 *    Description:  c functions which interface the python and c functions
 *
 *        Version:  0.2
 *        Created:  07/10/2014 14:54:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei-Yi Cheng (weiyi-bitw), ninpy.weiyi@gmail.com
 *   Organization:  Mount Sinai Medical Center
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>

/* ================= C functions ================= */

float log2f(float x) {
  return log(x)/log(2);
}

double log2d(double x) {
  return log(x)/log(2);
}


double mean(double *data, int numSamples) {
  int curSample;
  double mean = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    mean += data[curSample];
  }
  return mean / (double) numSamples;
}

double std(double *data, int numSamples) {
  int curSample;
  double m = mean(data, numSamples);
  double std = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    std += (data[curSample] - m)*(data[curSample] - m);
  }
  return sqrt(1/(double)(numSamples - 1) * std);
}

float meanf(float *data, int numSamples) {
  int curSample;
  float mean = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    mean += data[curSample];
  }
  return mean / (float) numSamples;
}

float stdf(float *data, int numSamples) {
  int curSample;
  float m = meanf(data, numSamples);
  float std = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    std += (data[curSample] - m)*(data[curSample] - m);
  }
  return sqrt(1/(float)(numSamples - 1) * std);
}

double meani(int *data, int numSamples) {
  int curSample;
  double mean = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    mean += (double)data[curSample];
  }
  return mean / (double) numSamples;
}

double stdi(int *data, int numSamples) {
  int curSample;
  double m = meani(data, numSamples);
  double std = 0;
	
  for (curSample = 0; curSample < numSamples; curSample++) {
    std += (double) (data[curSample] - m)*(data[curSample] - m);
  }
  return sqrt(1/(double)(numSamples - 1) * std);
}

/* Follows Daub et al, which contains mistakes; corrections based on spline descriptions on MathWorld pages */
double basisFunction(int i, int p, double t, const double *kVector, int numBins) {
  double d1, n1, d2, n2, e1, e2;
  if (p == 1) {
    if ((t >= kVector[i] && t < kVector[i+1] && 
	 kVector[i] < kVector[i+1]) ||
	(fabs(t - kVector[i+1]) < 1e-10 && (i+1 == numBins))) {
      return(1);
    }
    return(0);
  }
    
  d1 = kVector[i+p-1] - kVector[i];
  n1 = t - kVector[i];
  d2 = kVector[i+p] - kVector[i+1];
  n2 = kVector[i+p] - t;

  if (d1 < 1e-10 && d2 < 1e-10) {
    return(0);
  } else if (d1 < 1e-10) {
    e1 = 0;
    e2 = n2/d2*basisFunction(i+1, p-1, t, kVector, numBins);
  } else if (d2 < 1e-10) {
    e2 = 0;
    e1 = n1/d1*basisFunction(i, p-1, t, kVector, numBins);
  } else {
    e1 = n1/d1*basisFunction(i, p-1, t, kVector, numBins);
    e2 = n2/d2*basisFunction(i+1, p-1, t, kVector, numBins);
  }    
  
  /* sometimes, this value is < 0 (only just; rounding error); truncate */
  if (e1 + e2 < 0) {
    return(0);
  }
  return(e1 + e2);
}

void knotVector(double *v, int numBins, int splineOrder) {
  int nInternalPoints = numBins - splineOrder;
  int i;

  for (i = 0; i < splineOrder; ++i) {
    v[i] = 0;
  }
  for (i = splineOrder; i < splineOrder + nInternalPoints; ++i) {
    v[i] = (double)(i - splineOrder + 1)/(nInternalPoints + 1);
  }
  for (i = splineOrder + nInternalPoints; i < 2*splineOrder + nInternalPoints; ++i) {
    v[i] = 1;
  }
}

double max_d(const double *data, int numSamples) {
  int curSample;
  double curMax = data[0];
	
  for (curSample = 1; curSample < numSamples; curSample++) {
    if (data[curSample] > curMax) {
      curMax = data[curSample];
    }
  }
  return curMax;
}

double min_d(const double *data, int numSamples) {
  int curSample;
  double curMin = data[0];
	
  for (curSample = 1; curSample < numSamples; curSample++) {
    if (data[curSample] < curMin) {
      curMin = data[curSample];
    }
  }
  return curMin;
}

int maxi(const int *data, int numSamples) {
  int curSample;
  int curMax = data[0];
	
  for (curSample = 1; curSample < numSamples; curSample++) {
    if (data[curSample] > curMax) {
      curMax = data[curSample];
    }
  }
  return curMax;
}

int mini(const int *data, int numSamples) {
  int curSample;
  int curMin = data[0];
	
  for (curSample = 1; curSample < numSamples; curSample++) {
    if (data[curSample] < curMin) {
      curMin = data[curSample];
    }
  }
  return curMin;
}

void xToZ(const double *fromData, double *toData, int numSamples, int splineOrder, int numBins, double xMin, double xMax) {
  int curSample;
  
  if (xMin == -1 && xMax == -1) { /*then compute on the fly */
	  xMin = min_d(fromData, numSamples);
	  xMax = max_d(fromData, numSamples);
  } /*else use provided values */
  if(xMax == xMin) xMax = xMin + 1; // prevent "flat vector"
  for (curSample = 0; curSample < numSamples; curSample++) {
    /* toData[curSample] = (fromData[curSample] - xMin) * (numBins - splineOrder + 1) / (double) (xMax - xMin); */
    /* normalize to [0, 1] */
    toData[curSample] = (fromData[curSample] - xMin) / (double) (xMax - xMin);
  }
}

void findWeights(const double *x, const double *knots, double *weights, int numSamples, int splineOrder, int numBins, double rangeLeft, double rangeRight) {
  int curSample, curBin;
  double *z = (double*) calloc(numSamples, sizeof(double));

  xToZ(x, z, numSamples, splineOrder, numBins, rangeLeft, rangeRight);
  
  for (curSample = 0; curSample < numSamples; curSample++) {
    for (curBin = 0; curBin < numBins; curBin++) {
      weights[curBin * numSamples + curSample] = basisFunction(curBin, splineOrder, z[curSample], knots, numBins);
      //printf("%d|%f(%f)\t", curBin, weights[curBin * numSamples + curSample],z[curSample]);
    }
  }
  free(z);
}

void combineWeights(const double *wx, const double *wy, double *w, int numSamples, int numBins){
	int curSample, bx, by;
	for(bx = 0; bx < numBins; bx++){
	  for(by = 0; by < numBins; by++){
	    for(curSample = 0; curSample < numSamples; curSample++){
	    	w[ (bx*numBins + by) * numSamples + curSample] = wx[bx * numSamples + curSample] * wy[by * numSamples + curSample];
	    }
	  }
	}
}

double entropy1(const double *weights, int numSamples, int numBins) {
  int curBin, curSample;
  double H = 0, h;

  for (curBin = 0; curBin < numBins; curBin++) {
    h = 0;
    for(curSample = 0; curSample < numSamples; curSample++) h += weights[curBin * numSamples + curSample];
    h /= numSamples;
    if (h > 0) {
      H -= h * log2d(h);
    }
  }
  return H;
}

double entropy2(const double *wx, const double *wy, int numSamples, int numBins) {
  int curBinX, curBinY, curSample;
  double H = 0, h;

  for (curBinX = 0; curBinX < numBins; curBinX++) {
    for (curBinY = 0; curBinY < numBins; curBinY++) {
      h = 0;
      for(curSample = 0; curSample < numSamples; curSample++){
        h += wx[curBinX * numSamples + curSample] * wy[curBinY * numSamples + curSample];
      }
      h /= numSamples;
      if (h > 0) {
	H -= h * log2d(h);
      }
    }
  }
  return H;
}

double entropy2DiffBins(const double *wx, const double *wy, int numSamples, int binx, int biny){
	int curBinX, curBinY, curSample;
	double H = 0, h;

	for(curBinX = 0; curBinX < binx; curBinX++){
	  for(curBinY = 0; curBinY < biny; curBinY++){
		h = 0;
		for(curSample = 0; curSample < numSamples; curSample++){
			h += wx[curBinX*numSamples + curSample] * wy[curBinY * numSamples + curSample];
		}
		h /= numSamples;
		if(h > 0){
			H -= h * log2d(h);
		}
	  }
	}
	return H;
}

double entropy3(const double *wx, const double *wy, const double *wz, int numSamples, int numBins){
	int curBinX, curBinY, curBinZ, curSample;
	double H = 0, histVal = 0;
	double wxyz;
	for(curBinX = 0; curBinX < numBins; curBinX++){
	  for(curBinY = 0; curBinY < numBins; curBinY++){
	    for(curBinZ = 0; curBinZ < numBins; curBinZ++){
	      histVal = 0;
		for(curSample = 0; curSample < numSamples; curSample++){
		  wxyz = wx[curBinX*numSamples + curSample] * wy[curBinY*numSamples + curSample] * wz[curBinZ*numSamples + curSample];
		  histVal += wxyz;
		}
		histVal /= numSamples;
		if(histVal > 0){
			H -= histVal * log2d(histVal);
		}
	    }
	  }
	}
	return H;
}

double productMoment(const double *x, const double *y, int n){
	int i;
	double sumX=0, sumY=0, sumXY=0;
	for(i = 0; i < n; i++){
		sumX += x[i];
		sumY += y[i];
		sumXY += x[i] * y[i];
	}
	sumXY = sumXY * n - sumX * sumY;
	return sumXY;
}

//========================= export python function ===================================

double mi2(const double *x, const double *y, int n, int bin, int so, int norm, int negateMI){
  double *u = (double*) calloc(bin + so, sizeof(double));
  double *wx = (double*) calloc(bin * n, sizeof(double));
  double *wy = (double*) calloc(bin * n, sizeof(double));
  double e1x, e1y, mix, miy, largerMI, mi;

  knotVector(u, bin, so);
  findWeights(x, u, wx, n, so, bin, -1, -1);
  findWeights(y, u, wy, n, so, bin, -1, -1);
  e1x = entropy1(wx, n, bin);
  e1y = entropy1(wy, n, bin);
  mi = (e1x + e1y - entropy2(wx, wy, n, bin));

  if(norm == 1){
    mix = 2*e1x - entropy2(wx, wx, n, bin);
    miy = 2*e1y - entropy2(wy, wy, n, bin);
    largerMI = mix > miy ? mix:miy;
    if(largerMI == 0) largerMI = 1;
    mi /= largerMI;
  }
  if(negateMI == 1 && productMoment(x, y, n) < 0) mi = -mi;
  free(wx);
  free(wy);
  free(u);
  return mi;
}
/*
void mi2DiffBins(const double *x, const double *y, int *n, int *binx, int *biny, int *sox, int *soy, double *miOut, int *norm, int *negateMI){
	double *ux = (double*) calloc(*binx + *sox, sizeof(double));
	double *wx = (double*) calloc(*binx * *n, sizeof(double));
	double *uy = (double*) calloc(*biny + *soy, sizeof(double));
	double *wy = (double*) calloc(*biny * *n, sizeof(double));
	double e1x, e1y, mix, miy, largermi, mi;
	
	knotVector(ux, *binx, *sox);
	findWeights(x, ux, wx, *n, *sox, *binx, -1, -1);
	free(ux);
	knotVector(uy, *biny, *soy);
	findWeights(y, uy, wy, *n, *soy, *biny, -1, -1);
	free(uy);
	e1x = entropy1(wx, *n, *binx);
	e1y = entropy1(wy, *n, *biny);
	mi = (e1x + e1y - entropy2DiffBins(wx, wy, *n, *binx, *biny));
	
	if(*norm == 1){
		mix = 2*e1x - entropy2(wx, wx, *n, *binx);
		miy = 2*e1y - entropy2(wy, wy, *n, *biny);
		largermi = mix > miy? mix : miy;
		if(largermi == 0) largermi = 1;
		mi /= largermi;
	}
	if(*negateMI == 1 && productMoment(x, y, *n) < 0) mi = -mi;
	free(wx);
	free(wy);
	*miOut = mi;
}

void mi2vs1(const double *x, const double *y, const double *z, int *n, int *bin, int *so, double *mi, int *norm){
  double *u = (double*) calloc(*bin + *so, sizeof(double));
  double *wx = (double*) calloc(*bin * *n, sizeof(double));
  double *wy = (double*) calloc(*bin * *n, sizeof(double));
  double *wz = (double*) calloc(*bin * *n, sizeof(double));
  double *wxy = (double*) calloc(*bin * *bin * *n, sizeof(double));
  double e1xy, e1z, e2xyz, mixy, miz, largermi;

  knotVector(u, *bin, *so);
  findWeights(x, u, wx, *n, *so, *bin, -1, -1);
  findWeights(y, u, wy, *n, *so, *bin, -1, -1);
  findWeights(z, u, wz, *n, *so, *bin, -1, -1);
  combineWeights(wx, wy, wxy, *n, *bin);
  free(wx);
  free(wy);

  e1xy = entropy1(wxy, *n, *bin * *bin);
  e1z = entropy1(wz, *n, *bin);
  *mi = (e1xy + e1z - entropy2DiffBins(wxy, wz, *n, *bin * *bin, *bin));
  if(*norm == 1){
    mixy = 2*e1xy - entropy2(wxy, wxy, *n, *bin * *bin);
    miz = 2*e1z - entropy2(wz, wz, *n, *bin);
    largermi = mixy > miz ? mixy : miz;
    if(largermi == 0) largermi = 1;
    *mi /= largermi;
  }

  free(wxy);
  free(wz);
  free(u);
}
*/

void getAllMIWz(double *data, const double* vec, double *mi, int m, int n, int bin, int so, int norm, int negateMI){
  double *u = (double*) calloc(bin + so, sizeof(double));
  double *y = (double*) calloc(n, sizeof(double));
  double *wx = (double*) calloc(bin * n, sizeof(double));
  double *wy = (double*) calloc(bin * n, sizeof(double));
  int i, j;
  double e1x, e1y, mix, miy, largerMI;

  knotVector(u, bin, so);
  findWeights(vec, u, wx, n, so, bin, -1, -1);
  e1x = entropy1(wx, n, bin);
  mix = 2*e1x - entropy2(wx, wx, n, bin);

  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++) y[j] = data[j + i*n];
    //y = data[i];
    findWeights(y, u, wy, n, so, bin, -1, -1);
    e1y = entropy1(wy, n, bin);
    mi[i] = (e1x + e1y - entropy2(wx, wy, n, bin));
    if(norm == 1){
      largerMI = mix;
      miy = 2*e1y - entropy2(wy, wy, n, bin);
      if(miy > mix) largerMI = miy;
      if(largerMI == 0) largerMI = 1;
      mi[i] /= largerMI;
    }
    if(negateMI==1 && productMoment(y, vec, n) < 0) mi[i] = -mi[i];
  }

  free(wx);
  free(wy);
  free(u);
  free(y);
}
/*
void mi3(const double *x, const double *y, const double *z, int *n, int *bin, int *so, double *miOut){
  double *u = (double*) calloc(*bin + *so, sizeof(double));
  double *wx = (double*) calloc(*bin * *n, sizeof(double));
  double *wy = (double*) calloc(*bin * *n, sizeof(double));
  double *wz = (double*) calloc(*bin * *n, sizeof(double));
  double e1x, e1y, e1z, e2xy, e2yz, e2xz, mi;


  knotVector(u, *bin, *so);
  findWeights(x, u, wx, *n, *so, *bin, -1, -1);
  findWeights(y, u, wy, *n, *so, *bin, -1, -1);
  findWeights(z, u, wz, *n, *so, *bin, -1, -1);
  e1x = entropy1(wx, *n, *bin);
  e1y = entropy1(wy, *n, *bin);
  e1z = entropy1(wz, *n, *bin);
  e2xy = entropy2(wx, wy, *n, *bin);
  e2yz = entropy2(wy, wz, *n, *bin);
  e2xz = entropy2(wx, wz, *n, *bin);
  mi = entropy3(wx, wy, wz, *n, *bin) - e2xy - e2yz - e2xz + e1x + e1y + e1z;

  free(wx);
  free(wy);
  free(wz);
  free(u);
  *miOut = mi;
}

*/


/* =========== python interface ============== */

int not_doublematrix(PyArrayObject *mat){
    if(mat->descr->type_num != NPY_DOUBLE || mat->nd != 2){
        PyErr_SetString(PyExc_ValueError,
                "In not_doublematrix: array must be of type Float and 2 dimensional (n x m).");
        return 1;}
    return 0;
}


static PyObject *
basis_function(PyObject *self, PyObject* args){
    int i,p,numBins, listSize, j;
    double t;
    double *kVec;

    PyObject *listObj, *seq, *out;

    if(! PyArg_ParseTuple( args, "iidOi", &i, &p, &t, &listObj, &numBins )) return NULL;
    
    seq = PySequence_Fast(listObj, "Expected a sequence");
    listSize = PySequence_Size(listObj);
    kVec = (double*) calloc(listSize, sizeof(double));
    
    for(j = 0; j < listSize; j++)
        kVec[j] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, j));
    
    Py_DECREF(seq);
    out = Py_BuildValue("d", basisFunction(i, p, t, kVec, numBins));
    free(kVec);
    return out;
}
  
static PyObject*
find_weights(PyObject *self, PyObject* args){
    double *x, *knots, *weights;
    int numSamples, splineOrder=3, numBins=6, i;
    npy_intp w_dim[1] = {0};
    PyObject *xObj,  *seq;
    PyArrayObject *wObj;

    if(! PyArg_ParseTuple( args, "O|ii", &xObj, &numBins, &splineOrder )) return NULL;
    
    numSamples = PySequence_Size(xObj);
    x = (double*) calloc(numSamples, sizeof(double));
    knots = (double*) calloc(numBins + splineOrder, sizeof(double));
    //weights = (double*) calloc(numBins * numSamples, sizeof(double));

    seq = PySequence_Fast(xObj, "Expected a sequence");
    for(i = 0; i < numSamples; i++)
        x[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
    Py_DECREF(seq);

    knotVector(knots, numBins, splineOrder);
    
    w_dim[0] = numSamples * numBins;
    wObj = (PyArrayObject *)PyArray_SimpleNew(1, w_dim, PyArray_DOUBLE);
    weights = (double*) wObj->data;
    findWeights(x, knots, weights, numSamples, splineOrder, numBins, -1, -1);
      
    free(knots);
    free(x);
    
    return PyArray_Return(wObj);   
}

static PyObject*
entropy(PyObject *self, PyObject *args){
    double *x, *knots, *weights, H;
    int bins=6, so=3, i, n;
    PyObject *xObj, *seq, *out;

    if(! PyArg_ParseTuple( args, "O|ii", &xObj, &bins, &so )) return NULL;
    n = PySequence_Size(xObj);
    x = (double*) calloc(n, sizeof(double));
    knots = (double*) calloc(bins + so, sizeof(double));
    weights = (double*) calloc(bins * n, sizeof(double));
    
    seq = PySequence_Fast(xObj, "Expected a sequence");
    for(i = 0; i < n; i++)
        x[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq,i));
    Py_DECREF(seq);

    knotVector(knots, bins, so);
    findWeights(x, knots, weights, n, so, bins, -1, -1);
    H = entropy1(weights, n, bins);

    free(x);
    free(knots);
    free(weights);
    out = Py_BuildValue("d", H);
    return out;
}

static PyObject*
joint_entropy(PyObject *self, PyObject *args){
    double *x, *y, *knots, *wx, *wy, H;
    int bins=6, so=3, i, n;
    PyObject *xObj, *yObj, *seq, *seq2, *out;

    if(! PyArg_ParseTuple( args, "OO|ii", &xObj, &yObj, &bins, &so )) return NULL;

    n = PySequence_Size(xObj);
    if (n != PySequence_Size(yObj)){
        // something is wrong!
        return NULL;
    }
    x = (double*) calloc(n, sizeof(double));
    y = (double*) calloc(n, sizeof(double));
    knots = (double*) calloc(bins + so, sizeof(double));
    wx = (double*) calloc(bins * n, sizeof(double));
    wy = (double*) calloc(bins * n, sizeof(double));
    
    seq = PySequence_Fast(xObj, "Expected a sequence");
    seq2 = PySequence_Fast(yObj, "Expected a sequence");
    for(i = 0; i < n; i++){
        x[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq,i));
        y[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq2,i));
    }
    Py_DECREF(seq);
    Py_DECREF(seq2);

    knotVector(knots, bins, so);
    findWeights(x, knots, wx, n, so, bins, -1, -1);
    findWeights(y, knots, wy, n, so, bins, -1, -1);
    H = entropy2(wx, wy, n, bins);

    free(x);
    free(y);
    free(knots);
    free(wx);
    free(wy);
    out = Py_BuildValue("d", H);
    return out;
}

static PyObject*
mi(PyObject *self, PyObject *args){
    double *x, *y, MI;
    int bins=6, so=3, i, n, norm=1, negateMI=0;
    PyObject *xObj, *yObj, *seq, *seq2, *out;
    
    if(! PyArg_ParseTuple( args, "OO|iiii", &xObj, &yObj, &bins, &so, &norm, &negateMI )) return NULL;

    n = PySequence_Size(xObj);
    if (n != PySequence_Size(yObj)){
        // something is wrong!
        return NULL;
    }
    x = (double*) calloc(n, sizeof(double));
    y = (double*) calloc(n, sizeof(double));
    
    seq = PySequence_Fast(xObj, "Expected a sequence");
    seq2 = PySequence_Fast(yObj, "Expected a sequence");
    for(i = 0; i < n; i++){
        x[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq,i));
        y[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq2,i));
    }
    Py_DECREF(seq);
    Py_DECREF(seq2);

    MI = mi2(x, y, n, bins, so, norm, negateMI);

    free(x);
    free(y);

    out = Py_BuildValue("d", MI);
    return out;
}

static PyObject*
all_mi(PyObject *self, PyObject *args){
    double *data, *vec, *MI;
    int m, n, i, bins = 6, so = 3, norm = 1, negateMI = 1;
    npy_intp dim[1] = {0};
    PyArrayObject *dObj, *out;
    PyObject *vObj, *seq;

    if(! PyArg_ParseTuple( args, "OO|iiii", &dObj, &vObj, &bins, &so, &norm, &negateMI )) return NULL;
    
    if(not_doublematrix(dObj)) return NULL;

    m = dim[0] = dObj->dimensions[0];
    n = dObj->dimensions[1];

    if(n != PySequence_Size(vObj)) return NULL; // make sure input has compatible dimensions

    vec = (double*) calloc(n, sizeof(double));
    //MI = (double*) calloc(m, sizeof(double));

    seq = PySequence_Fast(vObj, "Expected a sequence");
    for(i = 0; i < n; i++)
        vec[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
    Py_DECREF(seq);
    
    data = (double*) dObj->data;
    out = (PyArrayObject*) PyArray_SimpleNew(1, dim, PyArray_DOUBLE);
    MI = (double*) out->data;
    
    getAllMIWz(data, vec, MI, m, n, bins, so, norm, negateMI);
 
    free(vec);
    return PyArray_Return(out);


    /*  
    numSamples = PySequence_Size(xObj);
    x = (double*) calloc(numSamples, sizeof(double));
    knots = (double*) calloc(numBins + splineOrder, sizeof(double));
    //weights = (double*) calloc(numBins * numSamples, sizeof(double));

    seq = PySequence_Fast(xObj, "Expected a sequence");
    for(i = 0; i < numSamples; i++)
        x[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
    Py_DECREF(seq);

    knotVector(knots, numBins, splineOrder);
    
    w_dim[0] = numSamples * numBins;
    wObj = (PyArrayObject *)PyArray_SimpleNew(1, w_dim, PyArray_DOUBLE);
    weights = (double*) wObj->data;
    findWeights(x, knots, weights, numSamples, splineOrder, numBins, -1, -1);
      
    free(knots);
    free(x);
    
    return PyArray_Return(wObj);   
    */
}


static PyMethodDef BSUtilMethods[] = 
{
    {"basis_function", basis_function, METH_VARARGS, "calculate basis function"},
    {"find_weights", find_weights, METH_VARARGS, "convert input vector into weight matrix"},
    {"entropy", entropy, METH_VARARGS, "calculate entropy given bins and spline order"},
    {"joint_entropy", joint_entropy, METH_VARARGS, "calculate joint entropy or two vectors given bins and spline order"},
    {"mi", mi, METH_VARARGS, "mutual information of two vector"},
    {"all_mi", all_mi, METH_VARARGS, "calculate mutual information between a vector and every row in a matrix"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
init_c_bsplinemi(void)
{
    (void) Py_InitModule("_c_bsplinemi", BSUtilMethods);
    import_array(); // This is EXTREMELY IMPORTANT!!!!
}




