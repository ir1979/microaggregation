#include <utility.h>

using namespace std;

class indexValueClass {
public:
	indexValueClass(){};
	double val;
	int idx;
};

struct indexValueCompareStruct {
	bool operator() (const indexValueClass &A, const indexValueClass &B) {
		return (A.val < B.val);
	}
} indexValueCompare;


//double** normalizeMatrix(double **data, double *&Ex, double *&stdDev, int NObjects, int NDims) {
//	Ex = calculateEx(data, NObjects, NDims); // new double[NDims];
//	stdDev = calculateStdDevN_1(data, NObjects, NDims);
//	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);
//
//	for (int j = 0; j < NDims; j++) 
//		if (stdDev[j]!=0) 
//			for (int i = 0; i < NObjects; i++)
//				allData[i][j] = (data[i][j] - Ex[j]) / stdDev[j];
//		else 
//			for (int i = 0; i < NObjects; i++)
//				allData[i][j] = (data[i][j] - Ex[j]);
//
//
//	//delete[] Ex;
//	//Ex = 0;
//	//delete[] stdDev;
//	//stdDev = 0;
//	return allData;
//}



precalculatedStats *precalculateAll(double** orig, int NRecords, int NDims) 
{
	precalculatedStats *PS = new precalculatedStats();
	int i, j, index;
	double t1, t2, t3;
	double* h = new double[1];

	int n_records_org = NRecords;
	int n_columns_org = NDims;

	PS->mean_org = new double[n_columns_org];

	PS->vec2_org = new double[n_columns_org];
	PS->vec4_org = new double[n_columns_org];

	PS->varm01_org = new double[n_columns_org];
	PS->varm2_org = new double[n_columns_org];

	PS->correlacions_org = new double[n_columns_org * n_columns_org];
	PS->mat11_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org*n_columns_org];
	PS->mat40_org = new double[n_columns_org * n_columns_org];
	PS->mat04_org = new double[n_columns_org * n_columns_org];
	PS->mat13_org = new double[n_columns_org * n_columns_org];
	PS->mat31_org = new double[n_columns_org * n_columns_org];
	PS->mat22_org = new double[n_columns_org * n_columns_org];
	PS->mat02_org = new double[n_columns_org * n_columns_org];
	PS->mat20_org = new double[n_columns_org * n_columns_org];


	PS->varm11_org = new double[n_columns_org * n_columns_org];
	PS->varr_org = new double[n_columns_org * n_columns_org];

	PS->histo_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];
	PS->varq_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];

	vector_moment_column(orig, n_records_org, n_columns_org, 1, PS->mean_org);
	vector_moment_central_column(orig, n_records_org, n_columns_org, 2, PS->vec2_org, PS->mean_org);
	vector_moment_central_column(orig, n_records_org, n_columns_org, 4, PS->vec4_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 1, PS->mat11_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 2, PS->mat22_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 4, 0, PS->mat40_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 4, PS->mat04_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 3, PS->mat13_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 3, 1, PS->mat31_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 0, PS->mat20_org, PS->mean_org);
	matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 2, PS->mat02_org, PS->mean_org);
	matrix_coef_correl_columns(orig, n_records_org, n_columns_org, PS->correlacions_org, PS->mean_org);

	//-------------------- Var(m01)
	calcular_varm01(PS->varm01_org, PS->vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m2)
	calcular_varm2(PS->varm2_org, PS->vec4_org, PS->vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m11)
	calcular_varm11(PS->varm11_org, PS->mat22_org, PS->mat11_org, n_columns_org, n_records_org);

	//-------------------- Var(r)

	for (i = 0; i < n_columns_org; i++) {
		for (j = 0; j < n_columns_org; j++) {
			index = i * n_columns_org + j;

			t1 = PS->mat22_org[index] / pow(PS->mat11_org[index], 2);

			t2 = PS->mat40_org[index] / pow(PS->mat20_org[index], 2) + PS->mat04_org[index] / pow(PS->mat02_org[index], 2);
			t2 += 2.0 * PS->mat22_org[index] / (PS->mat20_org[index] * PS->mat02_org[index]);
			t2 = t2 / 4.0;

			t3 = PS->mat31_org[index] / (PS->mat11_org[index] * PS->mat20_org[index]);
			t3 += PS->mat13_org[index] / (PS->mat11_org[index] * PS->mat02_org[index]);

			PS->varr_org[index] = pow(PS->correlacions_org[index], 2) / n_records_org * (t1 + t2 - t3);
		}
	}

	//-- Amount of appearances around each quantil --
	histograma_quantils(orig, n_columns_org, n_records_org, PS->histo_org, h);

	//-- var(Q) --
	calcular_varq(PS->varq_org, n_columns_org, n_records_org, PS->histo_org, h[0]);

	/*
	* --- Measures over masked data ---
	*/

	int enc = min(n_columns_org, 7);		//effective number of columns
	PS->nchoosekSaved = new int[8];
	for (int k = 1; k <= 7; k++) {
		PS->nchoosekSaved[k] = nchoosek(enc, k);
	}


	//producing all PS->combinations
	PS->comb1 = intAlloc2D_Fast(PS->nchoosekSaved[1],1);
	for (i = 0; i < PS->nchoosekSaved[1]; i++) {
		PS->comb1[i][0] = i;
	}

	PS->comb2 = intAlloc2D_Fast(PS->nchoosekSaved[2],2);

	int counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			PS->comb2[counter][0] = i1;
			PS->comb2[counter][1] = i2;
			counter++;
		}
	}

	PS->comb3 = intAlloc2D_Fast(PS->nchoosekSaved[3],3);

	counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			for (int i3 = i2 + 1; i3 < enc; i3++) {
				PS->comb3[counter][0] = i1;
				PS->comb3[counter][1] = i2;
				PS->comb3[counter][2] = i3;
				counter++;
			}
		}
	}

	PS->comb4 = intAlloc2D_Fast(PS->nchoosekSaved[4],4);

	counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			for (int i3 = i2 + 1; i3 < enc; i3++) {
				for (int i4 = i3 + 1; i4 < enc; i4++) {
					PS->comb4[counter][0] = i1;
					PS->comb4[counter][1] = i2;
					PS->comb4[counter][2] = i3;
					PS->comb4[counter][3] = i4;
					counter++;
				}
			}
		}
	}

	PS->comb5 = intAlloc2D_Fast(PS->nchoosekSaved[5],5);

	counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			for (int i3 = i2 + 1; i3 < enc; i3++) {
				for (int i4 = i3 + 1; i4 < enc; i4++) {
					for (int i5 = i4 + 1; i5 < enc; i5++) {
						PS->comb5[counter][0] = i1;
						PS->comb5[counter][1] = i2;
						PS->comb5[counter][2] = i3;
						PS->comb5[counter][3] = i4;
						PS->comb5[counter][4] = i5;
						counter++;
					}
				}
			}
		}
	}

	PS->comb6 = intAlloc2D_Fast(PS->nchoosekSaved[6],6);

	counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			for (int i3 = i2 + 1; i3 < enc; i3++) {
				for (int i4 = i3 + 1; i4 < enc; i4++) {
					for (int i5 = i4 + 1; i5 < enc; i5++) {
						for (int i6 = i5 + 1; i6 < enc; i6++) {
							PS->comb6[counter][0] = i1;
							PS->comb6[counter][1] = i2;
							PS->comb6[counter][2] = i3;
							PS->comb6[counter][3] = i4;
							PS->comb6[counter][4] = i5;
							PS->comb6[counter][5] = i6;
							counter++;
						}
					}
				}
			}
		}
	}

	PS->comb7 = intAlloc2D_Fast(PS->nchoosekSaved[7],7);

	counter = 0;
	for (int i1 = 0; i1 < enc; i1++) {
		for (int i2 = i1 + 1; i2 < enc; i2++) {
			for (int i3 = i2 + 1; i3 < enc; i3++) {
				for (int i4 = i3 + 1; i4 < enc; i4++) {
					for (int i5 = i4 + 1; i5 < enc; i5++) {
						for (int i6 = i5 + 1; i6 < enc; i6++) {
							for (int i7 = i6 + 1; i7 < enc; i7++) {
								PS->comb7[counter][0] = i1;
								PS->comb7[counter][1] = i2;
								PS->comb7[counter][2] = i3;
								PS->comb7[counter][3] = i4;
								PS->comb7[counter][4] = i5;
								PS->comb7[counter][5] = i6;
								PS->comb7[counter][6] = i7;
								counter++;
							}
						}
					}
				}
			}
		}
	}

	PS->range = new double[n_columns_org];
	double* rangeMin = new double[n_columns_org];
	double* rangeMax = new double[n_columns_org];

	for (int iCol = 0; iCol < n_columns_org; iCol++) {
		rangeMax[iCol] = rangeMin[iCol] = orig[0][iCol];
		for (int iRow = 1; iRow < n_records_org; iRow++) {
			if (orig[iRow][iCol] < rangeMin[iCol]) {
				rangeMin[iCol] = orig[iRow][iCol];
			} else if (orig[iRow][iCol] > rangeMax[iCol]) {
				rangeMax[iCol] = orig[iRow][iCol];
			}
		}
		double iRange = rangeMax[iCol] - rangeMin[iCol];
		PS->range[iCol] = iRange > 0 ? iRange : 1;
	}


	return PS;
}

int fact(int n) {
	if (n <= 1) {
		return 1;
	} else {
		return fact(n - 1) * n;
	}
}

int nchoosek(int n, int k) {
	if (n < k) {
		return 0;
	}

	return fact(n) / fact(k) / fact(n - k);
}

double** normalizeNS(double** dataOriginal, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataNormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataNormalized[i][j] = (dataOriginal[i][j] - mu1[j]) / (smu2[j] == 0 ? 1 : smu2[j]);
		}
	}

	return dataNormalized;
}

double** denormalizeNS(double** dataNormalized, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double* cmu1 = mean(dataNormalized,NRecords,NDims);
	double* csmu2 = stDev(dataNormalized,NRecords,NDims);

	double** dataDenormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataDenormalized[i][j] = (dataNormalized[i][j] - cmu1[j]) * smu2[j] / (csmu2[j] == 0 ? 1 : csmu2[j]) + mu1[j];
		}
	}

	return dataDenormalized;
}

double** denormalizeNS_Simple(double** dataNormalizedFull, double* mu1, double* smu2, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataDenormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			dataDenormalized[i][j] = (dataNormalizedFull[i][j]) * smu2[j] + mu1[j];
		}
	}

	return dataDenormalized;
}

int maxArray(int* arr, int n) {
	int M = arr[0];
	for (int i = 1; i < n; i++) {
		if (arr[i] > M) {
			M = arr[i];
		}
	}
	return M;
}

double** normalize_mM(double** data, long min, long max, int NRecords, int NDims) {
	int rowsCount = NRecords;
	int colsCount = NDims;

	double** dataNormalized = doubleAlloc2D_Fast(rowsCount,colsCount);

	double* minArray = minMatrix(data, rowsCount, colsCount);
	double* maxArray = maxMatrix(data, rowsCount, colsCount);

	for (int row = 0; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			dataNormalized[row][col] = round((data[row][col] - minArray[col]) / (maxArray[col] - minArray[col]) * (max - min)) + min;
		}
	}

	return dataNormalized;
}

double round(double number)
{
	return (number < 0.0)? ceil(number - 0.5) : floor(number + 0.5);
}


double* minMatrix(double** data, int rowsCount, int colsCount) {
	double* r = new double[colsCount];

	for (int j=0;j<colsCount;j++)
		r[j] = data[0][j];

	for (int row = 1; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			if (data[row][col] < r[col]) {
				r[col] = data[row][col];
			}
		}
	}

	return r;
}

double* maxMatrix(double** data, int rowsCount, int colsCount) {
	double* r = new double[colsCount];

	for (int j=0;j<colsCount;j++)
		r[j] = data[0][j];

	for (int row = 1; row < rowsCount; row++) {
		for (int col = 0; col < colsCount; col++) {
			if (data[row][col] > r[col]) {
				r[col] = data[row][col];
			}
		}
	}

	return r;
}


std::vector<std::vector<double> > doubleAlloc2D_vec(int r, int c, double initialValue)
{
	return std::vector<std::vector<double> >(r,std::vector<double>(c,initialValue));
}

double **doubleAlloc2D_rawFast(int r, int c)		// allocate n pts in dim
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c];		// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

double **doubleAlloc2D_Fast(int r, int c)		// allocate a 2dimensional array with r x dim zero  entries
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c]();		// allocate space for the array initialized to zero
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

double **doubleAlloc2D_Fast(int r, int c, double initialValue)		// allocate a 2dimensional array with r x dim zero  entries
{
	double ** pa = new double*[r];			// allocate pointers
	double  *  p = new double [r*c];		// allocate space for the array initialized to zero
	for (int i=0;i<r*c;i++)
		p[i] = initialValue;
	for (int i = 0; i < r; i++) 
		pa[i] = &(p[i*c]);

	return pa;
}

void doubleDelete2D_Fast(double **&pa)			// deallocate points
{
	delete [] pa[0];							// dealloc reserved storage
	delete [] pa;								// dealloc pointers
	pa = NULL;
}

std::vector<std::vector<int> > intAlloc2D_vec(int r, int c, int initialValue) {
	return std::vector<std::vector<int> > (r,std::vector<int>(c,initialValue));
}

int **intAlloc2D_rawFast(int r, int c) {

	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c];			// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

int **intAlloc2D_Fast(int r, int c) {
	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c]();		// allocate space for the array
	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;
}

int **intAlloc2D_Fast(int r, int c, int initialValue) {
	int ** pa = new int*[r];			// allocate pointers
	int  *  p = new int [r*c];		// allocate space for the array
	//std::memset (p,initialValue,r*c);
	for (int i=0;i<r*c;i++)
		p[i] = initialValue;

	for (int i = 0; i < r; i++) {
		pa[i] = &(p[i*c]);
	}
	return pa;

}

void intDelete2D_Fast(int **&pa) {
	delete [] pa[0];							// dealloc reserved storage
	delete [] pa;								// dealloc pointers
	pa = NULL;
}

int *intAlloc1D(int r) {
	return new int[r]();
}

int *intAlloc1D(int r, int initialValue) {
	int *t = new int[r];
	//std::memset(t,initialValue,r);
	std::fill(t, t + r, initialValue);

	return t;
}

double *doubleAlloc1D(int r) {
	return new double[r]();
}

double *doubleAlloc1D(int r, double initialValue) {
	double *t = new double[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}


float *floatAlloc1D(int r) {
	return new float[r]();
}

float *floatAlloc1D(int r, float initialValue) {
	float *t = new float[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}


bool *boolAlloc1D(int r) {
	return new bool[r]();
}

bool *boolAlloc1D(int r, bool initialValue) {
	bool *t = new bool[r];
	for (int i = 0; i < r; i++)
		t[i] = initialValue;

	return t;
}

void copy1D(int *dest, int *src, int n) {
	const int* s = src;
	int* d = dest;
	const int* const dend = dest + n;
	while ( d != dend )
		*d++ = *s++;	
	//std::memcpy(dest,src,n*sizeof(int));

	//for (int i = 0; i < n; i++)
	//	dest[i] = src[i];
}

void copy1D(double *dest, double *src, int n) {
	const double* s = src;
	double* d = dest;
	const double* const dend = dest + n;
	while (d != dend)
		*d++ = *s++;
	//std::memcpy(dest,src,n*sizeof(double));

	//for (int i = 0; i < n; i++)
	//	dest[i] = src[i];
}

void copy2D(double **dest, double **src, int r, int c) {
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			dest[i][j] = src[i][j];
}

void copy2DFrom1D(double **dest, double *src, int r, int c) {
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			dest[i][j] = src[i * c + j];
}

void print1D(int *info, int n, int offset) {
	for (int i = 0; i < n; i++)
		std::cout << info[i] + offset << " ";
}

void print1D(FILE *fp, int *info, int n, int offset) {
	for (int i = 0; i < n; i++)
		fprintf(fp, "%d ", info[i] + offset);
}

void print1D(double *info, int n, double offset) {
	for (int i = 0; i < n; i++)
		std::cout << info[i] + offset << " ";
}

void print2D(double **info, int r, int c, double offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			std::cout << info[i][j] + offset << " ";
		std::printf("\n");
	}
}

void print2D(FILE *fp, double **info, int r, int c, double offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			if (j<c-1)
				std::fprintf (fp, "%lf, ", info[i][j] + offset);
			else
				std::fprintf (fp, "%lf\n", info[i][j] + offset);

		//std::fprintf(fp, "\n");
	}
}


void print2D(int **info, int r, int c, int offset) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++)
			std::cout << info[i][j] + offset << " ";
		std::printf("\n");
	}
}

double* calculateNEx2(double** allData, int NObjects, int NDims) {
	double *NEx2 = new double[NDims];
	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < NObjects; i++)
			sum += pow2(allData[i][j]);
		NEx2[j] = sum;
	}
	return NEx2;
}

double* calculateEx(double **data, int NObjects, int NDims) {
	double *mean = new double[NDims];
	double revNObjects = 1.0 / NObjects;
	int i,j;
	for (j = 0; j < NDims; j++) {
		double sum = 0;
		for (i = 0; i < NObjects; i++)
			sum += data[i][j];
		mean[j] = sum * revNObjects;
	}
	return mean;
}

// consider rewriting
double* calculateExOnActiveData(double **data, int *active, int NObjects, int NDims) {
	double *mean = new double[NDims];

	int counter = 0;
	double sum;
	int i,j;


	sum = 0;
	for (i = 0; i < NObjects; i++)
		if (active[i]) {
			sum += data[i][0];
			counter++;
		}
		double revCounter = 1.0 / counter;

		mean[0] = sum * revCounter;


		for (j = 1; j < NDims; j++) {
			sum = 0;
			for (i = 0; i < NObjects; i++)
				if (active[i])
					sum += data[i][j];
			mean[j] = sum * revCounter;
		}

		return mean;
}

double calculateSST(double** allData, int NObjects, int NDims) {
	double *Ex = calculateEx(allData, NObjects, NDims); 
	double *NEx2 = calculateNEx2(allData, NObjects, NDims); 

	double SST = 0;
	for (int j = 0; j < NDims; j++) {
		SST += NEx2[j] - NObjects * Ex[j] * Ex[j];
	}

	delete[] Ex;
	delete[] NEx2;

	return SST;
}

double* calculateStdDevN_1(double** allData, int NObjects, int NDims) {
	double *Ex = calculateEx(allData, NObjects, NDims); // new double[NDims];
	double *stdDev = new double[NDims];

	for (int j = 0; j < NDims; j++) {
		double sum = 0;
		for (int i = 0; i < NObjects; i++)
			sum += (allData[i][j] - Ex[j])*(allData[i][j] - Ex[j]);
		stdDev[j] = sqrt(sum / (NObjects - 1));
	}

	delete[] Ex;
	Ex = 0;
	return stdDev;
}

double calculateDistance(double *x, double *y, int NDims) {
	double d = 0;
	for (int i = 0; i < NDims; i++)
		d += pow2C(x[i] - y[i]);
	return sqrt(d);
}

double calculateDistance2(double *x, double *y, int NDims) {
	double d = 0;
	for (int i = 0; i < NDims; i++)
		d += pow2C(x[i] - y[i]);
	return d;
}


double** normalizeMatrix(double **data, int NObjects, int NDims) {
	double *Ex = calculateEx(data, NObjects, NDims); // new double[NDims];
	double *stdDev = calculateStdDevN_1(data, NObjects, NDims);
	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);

	for (int j = 0; j < NDims; j++) 
		if (stdDev[j]!=0) 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - Ex[j]) / stdDev[j];
		else 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - Ex[j]);


	delete[] Ex;
	Ex = 0;
	delete[] stdDev;
	stdDev = 0;
	return allData;
}

double** normalizeMatrix(double **data, double *&mean, double *&stDev, int NObjects, int NDims) {
	mean = calculateEx(data, NObjects, NDims); // new double[NDims];
	stDev = calculateStdDevN_1(data, NObjects, NDims);
	double **allData = doubleAlloc2D_rawFast(NObjects, NDims);

	for (int j = 0; j < NDims; j++) 
		if (stDev[j]!=0) 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - mean[j]) / stDev[j];
		else 
			for (int i = 0; i < NObjects; i++)
				allData[i][j] = (data[i][j] - mean[j]);

	return allData;
}

double **denormalizeMatrix(double **centers, int *assignment, double *meanOrig, double *stdOrig, int NRecords, int NDims)
{
	double **result = doubleAlloc2D_Fast(NRecords,NDims);

	for (int i=0;i<NRecords;i++)
		copy1D(result[i],centers[assignment[i]],NDims);

	double *Ex = calculateEx(result, NRecords, NDims); // new double[NDims];
	double *stdDev = calculateStdDevN_1(result, NRecords, NDims);

	for (int j=0;j<NDims;j++)
		if (stdDev[j]!=0)
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])/stdDev[j]*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];
		else
			for (int i=0;i<NRecords;i++)
				result[i][j] = (result[i][j]-Ex[j])*((stdOrig[j]!=0)?stdOrig[j]:1)+meanOrig[j];

	return result;
}

int *getSortedIdxFast(double *A, int n) {

	//std::vector<indexValueClass> indexValue(n, indexValueClass()); 
	std::vector<indexValueClass> indexValue(n); 

	for (int idx = 0; idx < n; idx++) {
		indexValue[idx].val = A[idx];
		indexValue[idx].idx = idx;
	}

	std::sort (indexValue.begin(),indexValue.end(), indexValueCompare);

	int *result = intAlloc1D(n);

	for (int i=0;i<n;i++) {
		result[i] = indexValue[i].idx;
		A[i] = indexValue[i].val;
	}

	return result;

}



double sgn(double d) {
	if (d >= 0) {
		return 1;
	} else {
		return -1;
	}
}

void normalizeMm(double** org, double** msk, int n, int d) {
	double* minOrg = new double[d];
	double* maxOrg = new double[d];
	double* minMsk = new double[d];
	double* maxMsk = new double[d];
	for (int k = 0; k < d; k++) {
		minOrg[k] = maxOrg[k] = org[0][k];
		minMsk[k] = maxMsk[k] = msk[0][k];
	}

	for (int n1 = 1; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			if (minOrg[d1] > org[n1][d1]) {
				minOrg[d1] = org[n1][d1];
			}
			if (minMsk[d1] > msk[n1][d1]) {
				minMsk[d1] = msk[n1][d1];
			}
			if (maxOrg[d1] < org[n1][d1]) {
				maxOrg[d1] = org[n1][d1];
			}
			if (maxMsk[d1] < msk[n1][d1]) {
				maxMsk[d1] = msk[n1][d1];
			}
		}

	}

	for (int n1 = 0; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			org[n1][d1] = (org[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
			msk[n1][d1] = (msk[n1][d1] - minMsk[d1]) / (maxMsk[d1] - minMsk[d1]);
		}
	}
}

void normalizeMmGetMm(double** org, double** msk, int n, int d, double* maxOrg, double* minOrg) {
	//        double* minOrg = new double[d];
	//        double* maxOrg = new double[d];
	double* minMsk = new double[d];
	double* maxMsk = new double[d];
	for (int k = 0; k < d; k++) {
		minOrg[k] = maxOrg[k] = org[0][k];
		minMsk[k] = maxMsk[k] = msk[0][k];
	}

	for (int n1 = 1; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			if (minOrg[d1] > org[n1][d1]) {
				minOrg[d1] = org[n1][d1];
			}
			if (minMsk[d1] > msk[n1][d1]) {
				minMsk[d1] = msk[n1][d1];
			}
			if (maxOrg[d1] < org[n1][d1]) {
				maxOrg[d1] = org[n1][d1];
			}
			if (maxMsk[d1] < msk[n1][d1]) {
				maxMsk[d1] = msk[n1][d1];
			}
		}

	}

	for (int n1 = 0; n1 < n; n1++) {
		for (int d1 = 0; d1 < d; d1++) {
			org[n1][d1] = (org[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
			msk[n1][d1] = (msk[n1][d1] - minOrg[d1]) / (maxOrg[d1] - minOrg[d1]);
		}
	}
}

double** denormalizeMn(double** data01, int Nrecords, int NDims, double* maxMask, double* minMask) {
	double** data = doubleAlloc2D_Fast(Nrecords,NDims);

	for (int n1 = 0; n1 < Nrecords; n1++) {
		for (int d1 = 0; d1 < NDims; d1++) {
			data[n1][d1] = (maxMask[d1] - minMask[d1]) * data01[n1][d1] + minMask[d1];
		}
	}
	return data;
}

double** addArray(double** y, double** yp, int nrows, int ncols) {
	double**add = doubleAlloc2D_Fast(nrows,ncols);

	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			add[i][j] = y[i][j]+yp[i][j];
		}
	}
	return add;
}

double** addArray(double** y, double** yp, double c1, double c2, int nrows, int ncols) {

	double**add = doubleAlloc2D_Fast(nrows,ncols);

	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			add[i][j] = c1*y[i][j]+c2*yp[i][j];
		}
	}
	return add;
}

double** invertArrayPlus(double** y, int nrows, int ncols) {

	double **result =  doubleAlloc2D_Fast(nrows,ncols);
	for (int i=0;i<nrows;i++) {
		for (int j=0;j<ncols;j++) {
			result[i][j] = (y[i][j]);///(1+exp(-abs(y[i][j])));
		}
	}
	return result;
}


int knnsearch(double** data, int idx, int* exclude, int excludeCount, int nrows, int ncols) {
	double minDist = DOUBLE_MAX_VALUE;
	int nearIdx = -1;

	for (int i=0;i<nrows;i++) {
		bool flag = false;
		for (int cntr=0;cntr<excludeCount;cntr++) {
			if (i==exclude[cntr]) {
				flag=true;
				break;
			}
		}
		if (flag)
			continue;

		double tmpDist = dist_records(data, i, data, idx,ncols);

		if (tmpDist<minDist) {
			minDist = tmpDist;
			nearIdx = i;
		}
	}

	return nearIdx;
}


double** calculateCenters(double** allData, int* assignment, int NObjects, int NDims, int NClusters) {
	double** centers = doubleAlloc2D_Fast(NClusters,NDims);

	int* centerSizes = new int[NClusters];

	for (int i = 0; i < NClusters; i++) {
		centerSizes[i] = 0;
	}

	for (int i = 0; i < NObjects; i++) {
		centerSizes[assignment[i]]++;//=centerSizes[assignment[i]]+1;
	}
	for (int i = 0; i < NObjects; i++) {
		for (int j = 0; j < NDims; j++) {
			centers[assignment[i]][j] += allData[i][j];
		}
	}

	for (int i = 0; i < NClusters; i++) {
		for (int j = 0; j < NDims; j++) {
			centers[i][j] /= centerSizes[i];
		}
	}
	return centers;
}


void calculateDRIL(double** orig, double** mask, double &DR, double &IL, int NRecords, int NDims)
{	
	precalculatedStats *PS = precalculateAll(orig,NRecords,NDims);

	double* mu1Orig = calculateEx(orig,NRecords,NDims); 
	double* smu2Orig = calculateStdDevN_1(orig,NRecords,NDims);
	double **origNormalized = normalizeNS(orig,mu1Orig,smu2Orig,NRecords,NDims);

	double* mu1Mask = calculateEx(mask,NRecords,NDims);
	double* smu2Mask = calculateStdDevN_1(mask,NRecords,NDims);
	double **maskNormalized = normalizeNS(mask,mu1Mask,smu2Mask,NRecords,NDims);

	DR = calculateDisclosureRiskOptimized (orig, mask, origNormalized, maskNormalized, PS, NRecords, NDims); // last two args are for speedup
	IL = calculateInformationLossOptimized(orig, mask, PS, NRecords, NDims);

	delete []mu1Orig;
	delete []smu2Orig;

	delete []mu1Mask;
	delete []smu2Mask;

	doubleDelete2D_Fast(origNormalized);
	doubleDelete2D_Fast(maskNormalized);

	return;
}


double calculateInformationLossOptimized(double** orig, double** mask, precalculatedStats *PS, int NRecords, int NDims) 
{
	int n_columns_org, n_records_org;
	int n_columns_msk, n_records_msk;
	//        double* h = new double[1];

	double* mean_org = PS->mean_org;
	double* mean_msk;
	double* vec2_org = PS->vec2_org;
	double* vec2_msk;
	//        double* vec4_org=PS->vec4_org;
	double* varm01_org = PS->varm01_org;
	double* varm2_org = PS->varm2_org;

	double* correlacions_org = PS->correlacions_org;
	double* correlacions_msk;
	double* mat11_org = PS->mat11_org;
	double* mat11_msk;
	//        double* mat22_org=PS->mat22_org;
	//        double* mat40_org=PS->mat40_org;
	//        double* mat04_org=PS->mat04_org;
	//        double* mat13_org=PS->mat13_org;
	//        double* mat31_org=PS->mat31_org;
	//        double* mat02_org=PS->mat02_org;
	//        double* mat20_org=PS->mat20_org;
	double* varm11_org = PS->varm11_org;
	double* varr_org = PS->varr_org;
	//        double* histo_org=PS->histo_org;
	double* varq_org = PS->varq_org;

	//        int i, j, index;
	//        double t1, t2, t3;

	//        double* aux_msk_org;
	//        double* aux_v_org;
	//        double* aux_msk_msk;
	//        double* aux_v_msk;

	/*
	* Check for proper number of arguments
	*/

	n_records_org = NRecords;
	n_columns_org = NDims;
	//n_DATA_org = n_records_org * n_columns_org;

	n_records_msk = NRecords;
	n_columns_msk = NDims;
	//n_DATA_msk = n_records_msk * n_columns_msk;


	//mean_org = new double[n_columns_org];

	//vec2_org = new double[n_columns_org];
	//vec4_org = new double[n_columns_org];

	//varm01_org = new double[n_columns_org];
	//varm2_org = new double[n_columns_org];

	//correlacions_org = new double[n_columns_org * n_columns_org];
	//mat11_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org*n_columns_org];
	//mat40_org = new double[n_columns_org * n_columns_org];
	//mat04_org = new double[n_columns_org * n_columns_org];
	//mat13_org = new double[n_columns_org * n_columns_org];
	//mat31_org = new double[n_columns_org * n_columns_org];
	//mat22_org = new double[n_columns_org * n_columns_org];
	//mat02_org = new double[n_columns_org * n_columns_org];
	//mat20_org = new double[n_columns_org * n_columns_org];


	//varm11_org = new double[n_columns_org * n_columns_org];
	//varr_org = new double[n_columns_org * n_columns_org];

	//histo_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];
	//varq_org = new double[(100 / globalParameters::INC_QUANTIL + 1) * n_columns_org];

	//        vector_moment_column(orig, n_records_org, n_columns_org, 1, mean_org);
	//        vector_moment_central_column(orig, n_records_org, n_columns_org, 2, vec2_org, mean_org);
	//        vector_moment_central_column(orig, n_records_org, n_columns_org, 4, vec4_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 1, mat11_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 2, mat22_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 4, 0, mat40_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 4, mat04_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 1, 3, mat13_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 3, 1, mat31_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 2, 0, mat20_org, mean_org);
	//        matrix_moment_central_columns(orig, n_records_org, n_columns_org, 0, 2, mat02_org, mean_org);
	//        matrix_coef_correl_columns(orig, n_records_org, n_columns_org, correlacions_org, mean_org);

	//-------------------- Var(m01)
	//        calcular_varm01(varm01_org, vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m2)
	//        calcular_varm2(varm2_org, vec4_org, vec2_org, n_columns_org, n_records_org);

	//-------------------- Var(m11)
	//        calcular_varm11(varm11_org, mat22_org, mat11_org, n_columns_org, n_records_org);

	//-------------------- Var(r)

	//        for (i = 0; i < n_columns_org; i++) {
	//            for (j = 0; j < n_columns_org; j++) {
	//                index = i * n_columns_org + j;
	//
	//                t1 = mat22_org[index] / pow(mat11_org[index], 2);
	//
	//                t2 = mat40_org[index] / pow(mat20_org[index], 2) + mat04_org[index] / pow(mat02_org[index], 2);
	//                t2 += 2.0 * mat22_org[index] / (mat20_org[index] * mat02_org[index]);
	//                t2 = t2 / 4.0;
	//
	//                t3 = mat31_org[index] / (mat11_org[index] * mat20_org[index]);
	//                t3 += mat13_org[index] / (mat11_org[index] * mat02_org[index]);
	//
	//                varr_org[index] = pow(correlacions_org[index], 2) / n_records_org * (t1 + t2 - t3);
	//            }
	//        }

	//-- Amount of appearances around each quantil --
	//        histograma_quantils(orig, n_columns_org, n_records_org, histo_org, h);

	//-- var(Q) --
	//        calcular_varq(varq_org, n_columns_org, n_records_org, histo_org, h[0]);

	/*
	* --- Measures over masked data ---
	*/

	mean_msk = new double[n_columns_msk];
	vec2_msk = new double[n_columns_msk];
	mat11_msk = new double[n_columns_msk * n_columns_msk];
	correlacions_msk = new double[n_columns_msk * n_columns_msk];

	vector_moment_column(mask, n_records_msk, n_columns_msk, 1, mean_msk);
	vector_moment_central_column(mask, n_records_msk, n_columns_msk, 2, vec2_msk, mean_msk);
	matrix_moment_central_columns(mask, n_records_msk, n_columns_msk, 1, 1, mat11_msk, mean_msk);

	matrix_coef_correl_columns(mask, n_records_msk, n_columns_msk, correlacions_msk, mean_msk);

	//mexPrintf("\nProbabilistic comparison of files %s and %s (PIL(Q) PIL(m^0_1) PIL(m_2) PIL(m_{11}) PIL(r))\n",argv[1],argv[2]);

	//-- Quantiles information loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on quantiles PIL(Q): ");
	double PIL1 = loss_info_quantils(orig, mask, n_columns_org, n_records_org, varq_org);
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Means info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on means PIL(m^0_1): ");
	double PIL2 = loss_info_vector(mean_org, mean_msk, varm01_org, n_columns_org);
	//double PIL2 = 0;
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Variancies info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on variances PIL(m_2): ");
	double PIL3 = loss_info_vector(vec2_org, vec2_msk, varm2_org, n_columns_org);
	//double PIL3 = 0;
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Covariances info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on covariances PIL(m_{11}): ");
	double PIL4 = loss_info_mskatrix(mat11_org, mat11_msk, varm11_org, n_columns_org);
	//if(DETAILED_OUTPUT)mexPrintf("\n");
	//else mexPrintf(" ");

	//-- Correlation info. loss --
	//if(DETAILED_OUTPUT)mexPrintf("Average impact on Pearson's correlations PIL(r): ");
	double PIL5 = loss_info_mskatrix(correlacions_org, correlacions_msk, varr_org, n_columns_org);

	//if (DETAILED_OUTPUT) mexPrintf("\n");

	if (globalParameters::debugLevel>=10) {
		std::printf("\n PIL1: %f\tPIL2: %f\tPIL3: %f\tPIL4: %f\tPIL5: %f\t", PIL1, PIL2, PIL3, PIL4, PIL5);
	}
	//        std::printf (" IL: %f\t",(PIL1 + PIL2 + PIL3 + PIL4 + PIL5) / 5.0);
	return (PIL1 + PIL2 + PIL3 + PIL4 + PIL5) / 5.0;
	//if (DETAILED_OUTPUT) mexPrintf ("PIL: %lf\%", G_IL);
}


double moment_column(double** DATA, int n_records, int col, int n_columns, int r) {
	int i;
	double acum = 0.0;

	if (r == 1) {
		for (i = 0; i < n_records; i++) {
			acum += DATA[i][col];
		}
	} else {
		for (i = 0; i < n_records; i++) {
			acum += pow(DATA[i][col], r);
		}
	}
	acum /= n_records;
	return (acum);
}

//---------------------------------------------------
// Calcula r-essim moment central de la variable col
// Testejat amb Excel
//---------------------------------------------------
double moment_central_column(double** DATA, int n_records, int col, int n_columns, int r, double mean) {
	int i;
	//double mitja;
	double acum = 0.0;

	//mitja=moment_column(DATA,n_records,col,n_columns,1);

	for (i = 0; i < n_records; i++) {
		acum += pow(DATA[i][col] - mean, r);
	}
	acum /= n_records;
	return (acum);
}

//--------------------------------------------------------
// Calcula r1-r2 moment central de les variables col1, col2
// Testejat amb Excel
//--------------------------------------------------------
double moment_central_columns(double** DATA, int n_records, int col1, int col2, int n_columns, int r1, int r2, double* mean) {
	int i;
	/*
	* double mitja1, mitja2;
	*/
	double value1;
	double value2;
	double acum = 0.0;
	double resultat;

	//mitja1=moment_column(DATA,n_records,col1,n_columns,1);
	//mitja2=moment_column(DATA,n_records,col2,n_columns,1);

	for (i = 0; i < n_records; i++) {
		value1 = pow(DATA[i][col1] - mean[col1], r1);
		value2 = pow(DATA[i][col2] - mean[col2], r2);
		acum += (value1 * value2);
	}

	resultat = acum / n_records;
	return (resultat);
}

//-------------------------------------------------------------
// Calcula el coef. de correlacio de les variables col1, col2
// Testejat amb Excel
//-------------------------------------------------------------
double coef_correl_columns(double** DATA, int n_records, int col1, int col2, int n_columns, double* mean) {
	double mu11, mu20, mu02;
	double resultat;

	mu11 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 1, 1, mean);
	mu20 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 2, 0, mean);
	mu02 = moment_central_columns(DATA, n_records, col1, col2, n_columns, 0, 2, mean);

	resultat = mu11 / sqrt(mu20 * mu02);
	return (resultat);
}

//-------------------------------------------------------------
// Omple un vector amb la mitjana de cada column
//-------------------------------------------------------------
void vector_moment_column(double** DATA, int n_records, int n_columns, int r, double* mean) {
	int col;

	for (col = 0; col < n_columns; col++) {
		mean[col] = moment_column(DATA, n_records, col, n_columns, r);
	}
}

//-------------------------------------------------------------
// Omple un vector amb la varian\ca de cada column
//-------------------------------------------------------------
void vector_moment_central_column(double** DATA, int n_records, int n_columns, int r, double* variancies, double* mean) {
	int col;

	for (col = 0; col < n_columns; col++) {
		variancies[col] = moment_central_column(DATA, n_records, col, n_columns, r, mean[col]);
	}
}

//------------------------------------------------------------------------------
// Omple una matrix amb els r1-r2 moments centrals de cada parella de columns
//------------------------------------------------------------------------------
void matrix_moment_central_columns(double** DATA, int n_records, int n_columns, int r1, int r2, double *matrix, double* mean) {
	int col1, col2;

	for (col1 = 0; col1 < n_columns; col1++) {
		for (col2 = 0; col2 < n_columns; col2++) {
			matrix[col1 * n_columns + col2] = moment_central_columns(DATA, n_records, col1, col2, n_columns, r1, r2, mean);
		}
	}
}


//------------------------------------------------------------------------------
// Omple una matrix amb els coef. de correlacio de cada parella de columns
//------------------------------------------------------------------------------
void matrix_coef_correl_columns(double** DATA, int n_records, int n_columns, double* matrix, double* mean) {
	int col1, col2;

	for (col1 = 0; col1 < n_columns; col1++) {
		for (col2 = 0; col2 < n_columns; col2++) {
			matrix[col1 * n_columns + col2] = coef_correl_columns(DATA, n_records, col1, col2, n_columns, mean);
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm01(double* varm01, double* vec2, int dimensio, int n_records) {
	int i;

	for (i = 0; i < dimensio; i++) {
		varm01[i] = vec2[i] / (double) n_records;
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm2(double* varm2, double* vec4, double* vec2, int dimensio, int n_records) {
	int i;

	for (i = 0; i < dimensio; i++) {
		varm2[i] = (vec4[i] - pow(vec2[i], 2)) / (double) n_records;
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varm11(double* varm11, double* mat22, double* mat11, int dimensio, int n_records) {
	int i, j;
	int index;

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			index = i * dimensio + j;
			varm11[index] = (mat22[index] - pow(mat11[index], 2)) / (double) n_records;
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void histograma_quantils(double** DATA, int n_columns, int n_records, double* histograma, double* h) {
	double* vector;
	int i, j, k;
	double min, max;
	double quantil;
	int appearances;

	vector = new double[n_records];

	//--Per each variable
	for (i = 0; i < n_columns; i++) {
		crear_vector_orgrdenat(DATA, i, n_columns, n_records, vector);

		min = vector[0];
		max = vector[n_records - 1];

		h[0] = abs(max - min) / 1000.0;

		//-- Quantiles for each (globalParameters::INC_QUANTIL%) of the variable (I do not think the first or last)
		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {
			quantil = vector[(n_records - 1) * j / 100];

			appearances = 0;
			for (k = 0; k < n_records; k++) {
				if (vector[k] > (quantil - (h[0])) && vector[k] < (quantil + (h[0]))) {
					appearances++;
				} else if (appearances>0) {
                            break;
                            //k = k;
                }
			}

			histograma[n_columns * (j / globalParameters::INC_QUANTIL) + i] = (double) appearances;
			//	  mexPrintf("%d ",j/globalParameters::INC_QUANTIL);
			//      mexPrintf("%g %d\n",quantil, appearances);

		}
		//      mexPrintf("\n\n");

	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void calcular_varq(double* varq, int n_columns, int n_records, double* histo, double h) {
	int j, i;
	int index;

	for (i = 0; i < n_columns; i++) {
		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {
			index = n_columns * (j / globalParameters::INC_QUANTIL) + i;
			varq[index] = j * (100 - j) * n_records * pow(2 * h, 2) / (pow(histo[index], 2) * 100.0 * 100.0);
		}
	}
}

//-------------------------------------------------------------
//-------------------------------------------------------------
void crear_vector_orgrdenat(double** DATA, int col, int n_cols, int n_records, double* vector) {
	int i;

	/*
	* --- Primer copio la column sobre el vector ---
	*/
	for (i = 0; i < n_records; i++) {
		vector[i] = DATA[i][col];
	}

	/*
	* --- La ordeno ---
	*/
	
	std::sort(vector,vector+n_records);
	//quicksort(vector, 0, n_records - 1);
	
	
	/*
	* --- Verifico la ordenacio (aixo es podra treure) ---
	*/
	/*
	* for(i=0;i<(n_records-1);i++) if(vector[i]>vector[i+1]) {
	* mexPrintf("\nVector mal ordenat. Verificar quicksort!!!\n"); exit(0);
	* }
	*/
}

double p(double x) {
	double b1 = 0.319381530;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double pi = 3.14159265;

	double t, z, value;

	t = 1 / (1 + p * x);

	z = exp(-pow(x, 2) / 2.0) / sqrt(2.0 * pi);

	value = 0.5 - z * (b1 * t + b2 * pow(t, 2) + b3 * pow(t, 3) + b4 * pow(t, 4) + b5 * pow(t, 5));

	return (value);
}

//-�------------------------------------------------------
//--------------------------------------------------------
double loss_info_vector(double* v1, double* v2, double* var, int dimensio) {
	int i;
	double* perdues;
	double acum = 0.0;

	perdues = new double[dimensio];

	for (i = 0; i < dimensio; i++) {
		perdues[i] = 2 * 100.0 * p(abs(v1[i] - v2[i]) / sqrt(var[i]));
		acum += perdues[i];
	}

	return acum / (double) dimensio;
}

//-�------------------------------------------------------
//--------------------------------------------------------
double loss_info_mskatrix(double* mat1, double* mat2, double* var, int dimensio) {
	int i, j, index;
	double* perdues;
	double acum = 0.0;
	int n_DATA = 0;

	perdues = new double[dimensio * dimensio];

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			if (i != j) {
				index = i * dimensio + j;
				perdues[index] = 2 * 100.0 * p(abs(mat1[index] - mat2[index]) / sqrt(var[index]));

				acum += perdues[index];
				n_DATA++;
			}
		}
	}

	return acum / (double) n_DATA;
}

//-�------------------------------------------------------
//--------------------------------------------------------
void loss_info_mskatrix2(double* mat1, double* mat2, double* var, int dimensio) {
	int i, j, index;
	double* perdues;
	double acum = 0.0;
	int n_DATA = 0;

	perdues = new double[dimensio * dimensio];

	for (i = 0; i < dimensio; i++) {
		for (j = 0; j < dimensio; j++) {
			if (i != j) {
				index = i * dimensio + j;
				perdues[index] = 50.0 * abs(mat1[index] - mat2[index]);

				acum += perdues[index];
				n_DATA++;
			}
		}
	}
}

//-�------------------------------------------------------
//--------------------------------------------------------
double loss_info_quantils(double** DATA_org, double** DATA_msk, int n_columns, int n_records, double* varq_org) {
	double* vector_org;
	double* vector_msk;
	double dada_org, dada_msk, var_org;
	double prob;
	int i, j;
	double acum = 0.0;
	double count = 0;

	vector_org = new double[n_records];
	vector_msk = new double[n_records];

	for (i = 0; i < n_columns; i++) {
		crear_vector_orgrdenat(DATA_org, i, n_columns, n_records, vector_org);
		crear_vector_orgrdenat(DATA_msk, i, n_columns, n_records, vector_msk);

		for (j = globalParameters::INC_QUANTIL; j <= (100 - globalParameters::INC_QUANTIL); j += globalParameters::INC_QUANTIL) {

			dada_org = vector_org[(n_records - 1) * j / 100];
			dada_msk = vector_msk[(n_records - 1) * j / 100];
			var_org = varq_org[n_columns * (j / globalParameters::INC_QUANTIL) + i];

			prob = 2 * 100.0 * p(abs(dada_org - dada_msk) / sqrt(var_org));

			acum += prob;
			count++;
			//	  mexPrintf("%g %g %g %g -- ",dada_org,dada_msk,var_org,abs(dada_org-dada_msk));
			//	  mexPrintf("%g\n",prob);
		}
		//      mexPrintf("\n");
	}

	return acum / (double) count;
}

/*
* ---------------------------- quicksort -------------------------------
*/
void quicksort(double* vector, int inf, int sup) {
	int* k = new int[1];
	if (inf <= sup) {
		particio(vector, inf + 1, sup, vector[inf], k);
		intercanvi(vector, inf, k[0]);
		quicksort(vector, inf, k[0] - 1);
		quicksort(vector, k[0] + 1, sup);
	}
}
/*
* ----------------------------- intercanvi -----------------------------
*/

void intercanvi(double* arr, int a, int b) {
	double temp;
	temp = arr[a];
	arr[a] = arr[b];
	arr[b] = temp;
}
/*
* ----------------------------- particio --------------------------------
*/

void particio(double* vector, int inf, int sup, double x, int* k) {
	int k2;

	k[0] = inf - 1;
	k2 = sup + 1;
	while (k2 != (k[0] + 1)) {
		if (vector[k[0] + 1] <= x) {
			(k[0])++;
		} else if (vector[k2 - 1] >= x) {
			k2--;
		} else {
			intercanvi(vector, k[0] + 1, k2 - 1);
			(k[0])++;
			k2--;
		}
	}
}

/*
* ------------------------- cerca dicotomica -----------------------
*/
void cerca_dicotomica(double* vector, int mida_vector, double valor_buscat, int* pos) {
	int inf, sup, mig;

	inf = 0;
	sup = mida_vector - 1;

	while (inf != (sup + 1)) {
		mig = (sup + inf) / 2;

		if (vector[mig] <= valor_buscat) {
			inf = mig + 1;
		} else {
			sup = mig - 1;
		}
	}
	pos[0] = sup;
}


double calculateDisclosureRiskOptimized(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int NRecords, int NDims) {

	int n_data, n_columns, n_regs;
	//int i,j;
	double* dlds = new double[7];
	double* ids = new double[10];
	double dld, id;

	n_regs = NRecords;
	n_columns = NDims;
	n_data = n_regs * n_columns;

	int* match_all = new int[n_regs];

	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
	for (int kk = 0; kk < n_regs; kk++) {
		match_all[kk] = kk;
	}



	if (globalParameters::ignoreDLD==false) {
		calculate_dld_combination_optimized(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS);//comb1, comb2, comb3, comb4, comb5, comb6, comb7);
		dld = 0;
		for (int j = 0; j < min(7, n_columns); j++) {
			//std::printf ("%f ", dlds[j]); ///delme
			dld += dlds[j];
		}
	}
	else {
		dld = 0;
	}


	if (globalParameters::standard01) {
		//Normalize using max and min

		double **orig01 = doubleAlloc2D_Fast(n_regs, n_columns);
		double **mask01 = doubleAlloc2D_Fast(n_regs, n_columns);

		copy2D(orig01,orig,n_regs,n_columns);
		copy2D(mask01,mask,n_regs,n_columns);

		normalizeMm(orig01, mask01, n_regs, n_columns);

		//calculate id on 0-1 normalized data
		calculate_id(orig01, mask01, n_data, n_columns, match_all, ids);

		doubleDelete2D_Fast (orig01);
		doubleDelete2D_Fast (mask01);
	} else if (globalParameters::sdid) {
		calculate_sdid_relative(orig, mask, n_data, n_columns, match_all, ids);
	} else
	{
		calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
	}

	//        calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);

	id = 0;
	for (int j = 0; j < 10; j++) {
		id += ids[j];
	}
	id /= 10.0;

	///////////////////

	if (globalParameters::debugLevel>=10) {
		std::printf ("\tdld: %f ", dld);
		std::printf ("\tid: %f\n", id);
	}

	delete []dlds;
	delete []ids;
	delete []match_all;

	double dRisk;

	if (globalParameters::ignoreDLD)
		dRisk = id;
	else
		dRisk = (dld + id) / 2.0;

	return dRisk;

}



double calculateDisclosureRiskOptimizedRelative(double** orig, double** mask, double** origNormal, double** maskNormal, precalculatedStats *PS, int nrows, int ncols) {

	int n_data, n_columns, n_regs;
	//int i,j;
	double* dlds = new double[7];
	double* ids = new double[10];
	double dld, id;

	n_regs = nrows;
	n_columns = ncols;
	n_data = n_regs * n_columns;

	int* match_all = new int[n_regs];

	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
	for (int kk = 0; kk < n_regs; kk++) {
		match_all[kk] = kk;
	}

	//do_matching(data,n_regs,n_columns,mask,n_regs,n_columns,match_all,n_columns);/*-- last parameter matching n_vbles --*/
	//if(DEBUG)printf("matching done\n");

	//producing all combinations
	//int enc = min(n_columns, 7);		//effective number of columns
	//int** comb1 = PS->comb1;

	//int** comb2 = PS->comb2;

	//int** comb3 = PS->comb3;

	//int** comb4 = PS->comb4;

	//int** comb5 = PS->comb5;

	//int** comb6 = PS->comb6;

	//int** comb7 = PS->comb7;

	//-- with standardized data
	//        for (int ii=0;ii<3;ii++)
	//        {
	//            for (int jj=0;jj<5;jj++)
	//                std::printf("%f, ", origNormal[ii][jj]);
	//            System.out.println();
	//        }
	//
	//        System.out.println();
	//        System.out.println();
	//
	//        for (int ii=0;ii<3;ii++)
	//        {
	//            for (int jj=0;jj<5;jj++)
	//                std::printf("%f, ", maskNormal[ii][jj]);
	//            System.out.println();
	//        }


	calculate_dld_combination_optimized(origNormal, maskNormal, n_data, n_columns, match_all, dlds, PS);//comb1, comb2, comb3, comb4, comb5, comb6, comb7);
	//calculate_dld(data,mask,n_data,n_columns,match_all,dlds);

	/*
	* //-- I reset the data without normalize for(j=0;j<n_data;j++)
	* data[j]=orig[j]; for(j=0;j<n_data;j++) mask[j]=mask_copy[j];
	*
	* //Normalize using max and min
	* normalizeMm(data,mask,n_regs,n_columns);
	* calculate_id(data,mask,n_data,n_columns,match_all,ids);
	* if(DEBUG)printf("id calculated\n");
	*/

	//        for (int ii=0;ii<20;ii++)
	//        {
	//            for (int jj=0;jj<2;jj++)
	//                std::printf("%f ", origNormal[ii][jj]);
	//            std::printf ("\n");
	//        }
	//
	//        std::printf ("\n");
	//        std::printf ("\n");
	//                
	//        for (int ii=0;ii<20;ii++)
	//        {
	//            for (int jj=0;jj<2;jj++)
	//                std::printf("%f ", maskNormal[ii][jj]);
	//            std::printf ("\n");
	//        }

	if (globalParameters::standard01) {
		//Normalize using max and min
		double** orig01 = cloneArray(orig,n_regs, n_columns);
		double** mask01 = cloneArray(mask,n_regs, n_columns);

		normalizeMm(orig01, mask01, n_regs, n_columns);


		//calculate id on 0-1 normalized data
		calculate_id(orig01, mask01, n_data, n_columns, match_all, ids);
	} else {
		calculate_id_relative(orig, mask, n_data, n_columns, match_all, ids);
		//            calculate_id_relative(origNormal, maskNormal, n_data, n_columns, match_all, ids);
	}


	dld = 0;
	for (int j = 0; j < min(7, n_columns); j++) {
		dld += dlds[j];
	}

	id = 0;
	for (int j = 0; j < 10; j++) {
		id += ids[j];
	}
	id /= 10.0;

	///////////////////

	if (globalParameters::debugLevel>=10) {
		std::printf("dld: %f\t", dld);
		std::printf("id: %f\t", id);
	}

	return (dld + id) / 2.0;

}

void do_matching_combination(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}

	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = 999999;
		quin = -1;

		for (j = 0; j < n_regs2; j++) {
			distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			//printf ("d(%d,%d)=%f ",i,j,distancia);
			if (quin == -1) {
				min_distancia = distancia;
				quin = j;
			} else {
				if (distancia < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
					min_distancia = distancia;
					quin = j;
				}
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -i;
		}
	}
}

void do_matching_combination_optimized(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}

	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = 999999;
		quin = -1;

		for (j = 0; j < n_regs2; j++) {
			//distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			if (quin == -1) {
				min_distancia = distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
				quin = j;
			} else {
				//if ((distancia = dist_records_combination_with_threshold(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles, min_distancia)) < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
				if ((distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles)) < min_distancia) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
					min_distancia = distancia;
					quin = j;
				}
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -i;
		}
	}
}

void do_matching_combination_optimized2(double** data1, int n_regs1, int n_columns1, double** data2, int n_regs2, int n_columns2, int *resultat, int n_vbles, int* vbles) {
	int i, j;
	//double EPSILON = 1e-15;
	double distancia, min_distancia;
	int quin;

	//  printf("%d\n%d\n%d\n%d\n",n_regs1,n_columns1,n_regs2,n_columns2);
	//  printf("%d\n",n_vbles);

	/*
	* -- Inicialitzo --
	*/
	for (i = 0; i < n_regs1; i++) {
		resultat[i] = -1;
	}



	/*
	* -- For all registration date2 --
	*/
	for (i = 0; i < n_regs1; i++) {
		min_distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, i, n_vbles, vbles);
		quin = i;

		for (j = 0; j < n_regs2; j++) {
			if (j == i) {
				continue;
			}
			distancia = dist_records_combination(data1, n_columns1, data2, n_columns2, i, j, n_vbles, vbles);
			if ((distancia < min_distancia) || (distancia == min_distancia && j < i)) {// && ((min_distancia-distancia)>EPSILON || j==i)) {     //This is due to numerical hardware problems in very low numbers
				quin = -1;
				break;
			}
		}

		if (quin != -1) {
			//                if (i==0)
			//                    std::printf("%d\t%d\t%g\n",quin,i,min_distancia);
			resultat[i] = quin;
		} else {
			resultat[i] = -1;
		}
	}
}

/*----------------------------------------------------------------------------
Returns the Euclidean distance squared between data1 and register of i1
the register of i2 data2
----------------------------------------------------------------------------*/
double dist_records(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles) {
	double distancia = 0;
	int i;

	for (i = 0; i < n_vbles; i++) {
		distancia += (data1[i1][i] - data2[i2][i]) * (data1[i1][i] - data2[i2][i]);
	}
	return (distancia);
}

double dist_records(double** data1, int i1, double** data2, int i2,int ncols) {
	double distancia = 0;
	int i;

	for (i = 0; i < ncols; i++) {
		distancia += (data1[i1][i] - data2[i2][i]) * (data1[i1][i] - data2[i2][i]);
	}
	return sqrt(distancia);
}


/*
* ----------------------------------------------------------------------------
* Returns the Euclidean distance squared between data1 and register of i1
* the register of i2 data2, the combination is stores in vbles array
----------------------------------------------------------------------------
*/
double dist_records_combination(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles, int* vbles) {
	double distancia = 0;
	int i;
	double tmpValue;

	for (i = 0; i < n_vbles; i++) {
		tmpValue = (data1[i1][vbles[i]] - data2[i2][vbles[i]]);
		distancia += tmpValue * tmpValue;
	}
	return distancia;
}

double dist_records_combination_with_threshold(double** data1, int n_columns1, double** data2, int n_columns2, int i1, int i2, int n_vbles, int* vbles, double threshold) {
	double distancia = 0;
	int i;
	double tmpValue;

	for (i = 0; i < n_vbles; i++) {
		tmpValue = (data1[i1][vbles[i]] - data2[i2][vbles[i]]);
		distancia += tmpValue * tmpValue;
		if (distancia >= threshold) {
			return DOUBLE_MAX_VALUE;
		}
	}
	return distancia;
}

void calculate_dld_combination_optimized(double** data, double** data_masked, int n_data, int n_columns, int *match_all, double *dld, precalculatedStats *PS) //int** comb1, int** comb2, int** comb3, int** comb4, int** comb5, int** comb6, int** comb7) {
{
	int n_regs;
	int* match;
	int i, j, k, equal;
	double percentatge;
	int counter, finalCounter;

	n_regs = n_data / n_columns;

	match = new int[n_regs];

	counter = 0;
	finalCounter = (1 << min(7, n_columns)) - 1;		//-1 because no attribute is not allowed for DLD

	for (i = 1; i <= min(7, n_columns); i++) {
		equal = 0;
		for (k = 0; k < PS->nchoosekSaved[i]; k++, counter++) {
			switch (i) {
			case 1:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb1[k]);
				break;
			case 2:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb2[k]);
				break;
			case 3:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb3[k]);
				break;
			case 4:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb4[k]);
				break;
			case 5:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb5[k]);
				break;
			case 6:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb6[k]);
				break;
			case 7:
				do_matching_combination_optimized2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match, i, PS->comb7[k]);
				break;
			default:
				break;
			}

			for (j = 0; j < n_regs; j++) {
				if (match[j] == j) {
					equal++;
				}
			}
		}

		//std::printf ("eq=%d\n\n",equal);
		percentatge = (double) equal * 100.0 / (double) n_regs / finalCounter;
		dld[i - 1] = percentatge;
	}
}

/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
void calculate_id(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	i = 0;
	p = 1.0;
	while (i < 10) {
		//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
		percentatge = look_at_intervals2(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p);

		ids[i] = percentatge;
		//            std::printf("%d -> %f\t", i, percentatge);
		i++;
		p += 1.0;
	}
}


void calculate_sdid_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	double *sd = calculateStdDevN_1(data_masked,n_regs,n_columns);

	if (globalParameters::sdidRatio<=0) {
		i = 0;
		p = 1.0;
		while (i < 10) {
			//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
			percentatge = look_at_intervals2_SDIDrelative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p,sd);

			ids[i] = percentatge;
			i++;
			p += 1.0;
		}
	}
	else {
		percentatge = look_at_intervals2_SDIDrelative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, globalParameters::sdidRatio,sd);
		for (int i=0;i<10;i++)
			ids[i] = percentatge;
	}

	delete []sd;
}


/*-----------------------------------------------------------------------------
-----------------------------------------------------------------------------*/
void calculate_id_relative(double** data, double** data_masked, int n_data, int n_columns, int* match_all, double* ids) {
	double p;
	double percentatge;
	int i, n_regs;

	n_regs = n_data / n_columns;

	i = 0;
	p = 1.0;
	while (i < 10) {
		//percentatge=look_at_intervals(data,n_regs,n_columns,data_masked,n_regs,n_columns,match_all,p);
		percentatge = look_at_intervals2_relative(data, n_regs, n_columns, data_masked, n_regs, n_columns, match_all, p);

		ids[i] = percentatge;
		i++;
		p += 1.0;
	}

}

double look_at_intervals2(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;


	for (vble = 0; vble < n_vbles1; vble++) 
	{
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];
			v_min = valor_masked - p * valor_masked / 100.0;
			v_max = valor_masked + p * valor_masked / 100.0;
			//                v_min = valor_masked - p / 100.0;
			//                v_max = valor_masked + p / 100.0;

			if (valor_original >= v_min && valor_original <= v_max) {
				within++;
			} else {
				out++;
			}
		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}

double look_at_intervals2_SDIDrelative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p, double *SDmask) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;

	for (vble = 0; vble < n_vbles1; vble++)  {
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];

			////////////////////////////////////////////////

			v_min = valor_masked - p * abs(SDmask[vble] / 100.0);
			v_max = valor_masked + p * abs(SDmask[vble] / 100.0);
			//                v_min = valor_masked - p / 100.0;
			//                v_max = valor_masked + p / 100.0;

			if (valor_original >= v_min && valor_original <= v_max) {
				within++;
			} else {
				out++;
			}

			///////////////////////////////////////////////////

		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}


double look_at_intervals2_relative(double** data1, int n_regs1, int n_vbles1, double** data2, int n_regs2, int n_vbles2, int *match, double p) {
	int vble, i, i_ori;
	double valor_masked, valor_original, v_min, v_max;
	int pos;
	int within = 0, out = 0;

	for (vble = 0; vble < n_vbles1; vble++)  {
		for (i = 0; i < n_regs2; i++) {
			valor_original = data1[i][vble];
			valor_masked = data2[i][vble];

			////////////////////////////////////////////////

			v_min = valor_original - p * abs(valor_original / 100.0);
			v_max = valor_original + p * abs(valor_original / 100.0);
			//                v_min = valor_masked - p / 100.0;
			//                v_max = valor_masked + p / 100.0;

			if (valor_masked >= v_min && valor_masked <= v_max) {
				within++;
			} else {
				out++;
			}

			///////////////////////////////////////////////////

		}
	}

	return (100.0 * (double) within / ((double) within + (double) out));
}


double* mean(double** data, int rowsCount, int colsCount) 
{

	double* result = new double[colsCount];

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			result[j] += data[i][j];
		}
	}

	for (int j = 0; j < colsCount; j++) {
		result[j] /= rowsCount;
	}

	return result;
}

double* stDev(double** data, int rowsCount, int colsCount) 
{

	double* result = new double[colsCount];
	double* result2 = new double[colsCount];

	for (int i = 0; i < rowsCount; i++) {
		for (int j = 0; j < colsCount; j++) {
			result[j] += data[i][j];
			result2[j] += data[i][j] * data[i][j];
		}
	}

	for (int j = 0; j < colsCount; j++) {
		result[j] = sqrt((result2[j] - result[j] * result[j] / rowsCount) / (rowsCount - 1));
	}

	return result;
}

double **cloneArray(double **data, int nrows, int ncols)
{
	double **result = doubleAlloc2D_Fast(nrows,ncols);
	copy2D(result,data,nrows,ncols);
	return result;
}



//int *getSortedIdxFast(double *A, int n) {
//
//	//std::vector<indexValueClass> indexValue(n, indexValueClass()); 
//	std::vector<indexValueClass> indexValue(n); 
//
//	for (int idx = 0; idx < n; idx++) {
//		indexValue[idx].val = A[idx];
//		indexValue[idx].idx = idx;
//	}
//
//	std::sort (indexValue.begin(),indexValue.end(), indexValueCompare);
//
//	int *result = intAlloc1D(n);
//
//	for (int i=0;i<n;i++) {
//		result[i] = indexValue[i].idx;
//		A[i] = indexValue[i].val;
//	}
//
//	return result;
//
//}
//
//
//


//inline double calculatedistance(double *x, double *y, int ndims) {
//	double d = 0;
//	for (int i = 0; i < ndims; i++)
//		d += pow2c(x[i] - y[i]);
//	return sqrt(d);
//}

//inline double calculateDistance2(double *x, double *y, int NDims) {
//	double d = 0;
//	for (int i = 0; i < NDims; i++)
//		d += pow2C(x[i] - y[i]);
//	return d;
//}

