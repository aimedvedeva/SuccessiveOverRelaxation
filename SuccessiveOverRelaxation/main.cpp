#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

using namespace std;

ifstream fileOut;


void GetMatrix(int numOfVar, int numOfEq, double *matrix/*out*/) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			fileOut >> matrix[i * numOfVar + j];
		}
	}
}

void GetMartixSize(int *numOfVar/*out*/, int *numOfEq/*out*/) {
	fileOut >> *numOfEq >> *numOfVar;
	return;
}

void initMatrixByZero(double *matrix, int numOfVar, int numOfEq) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			matrix[i * numOfVar + j] = 0;
		}
	}
	return;
}

void InitDiagMatrixFromCoefMatrix(double *matrix, double *D, int numOfVar, int numOfEq) {
	initMatrixByZero(D, numOfVar, numOfEq);

	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			if (i == j) {
				D[i * numOfVar + j] = matrix[i * numOfVar + j];
			}
		}
	}
	return;
}

void InitUpperTriagMatrWithoutDiag(double *matrix, double *R, int numOfVar, int numOfEq) {
	initMatrixByZero(R, numOfVar, numOfEq);
	
	for (int i = 0; i < numOfEq - 1; i++) {
		for (int j = i + 1; j < numOfVar; j++) {
			R[i * numOfVar + j] = matrix[i * numOfVar + j];
		}
	}
	return;
}

void InitLowerTriagMatrWithoutDiag(double *matrix, double *L, int numOfVar, int numOfEq) {
	initMatrixByZero(L, numOfVar, numOfEq);

	for (int i = 1; i < numOfEq; i++) {
		for (int j = 0; j < i; j++) {
			L[i * numOfVar + j] = matrix[i * numOfVar + j];
		}
	}
	return;
}

void MultMatrixByNum(double *matrix, int numOfColumns, int numOfRows, double num) {
	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrix[i * numOfColumns + j] *= num;
		}
	}
}

double *MultMatrixAB(double *matrixA, int numOfVariablesA, int numOfEquationsA, double *matrixB, int numOfVariablesB, int numOfEquationsB) {
	if (numOfVariablesA != numOfEquationsB) {
		return NULL;
	}

	double *newMatrix = new double[numOfEquationsA * numOfVariablesB];

	for (int i = 0; i < numOfEquationsA; i++) {
		for (int j = 0; j < numOfVariablesB; j++) {
			double sum = 0;
			for (int k = 0; k < numOfVariablesA; k++) {
				sum += matrixA[i * numOfVariablesA + k] * matrixB[k * numOfVariablesB + j];
			}
			newMatrix[i * numOfVariablesB + j] = sum;
		}
	}

	return newMatrix;
}

void InitInvertibleDiagMatrix(double *D, int numOfVar, int numOfEq) {
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			if (i == j) {
				D[i * numOfVar + j] = 1 / D[i * numOfVar + j];
			}
		}
	}
	return;
}

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
	}
	sum = sqrt(sum);
	return sum;
}

void SubtractMatrix(double *matrixA, double *matrixB, int numOfColumns, int numOfRows) {

	for (int i = 0; i < numOfRows; i++) {
		for (int j = 0; j < numOfColumns; j++) {
			matrixA[i * numOfColumns + j] -= matrixB[i * numOfColumns + j];
		}
	}
}

void CopyMatrix(double *to, double *from, int numOfRows, int numOfColumns, int fromColumn, int fromRow) {
	for (int i = fromRow; i < numOfRows; i++) {
		for (int j = fromColumn; j < numOfColumns; j++) {
			to[i * numOfColumns + j] = from[i * numOfColumns + j];
		}
	}
}

bool isEnoughIteration(double *prevGuess, double *guess, int lenght, double approximation, double normC) {
	double *result = new double[lenght];

	CopyMatrix(result, guess, lenght, 1, 0, 0);
	SubtractMatrix(result, prevGuess, 1, lenght);

	double norm = EuclideanNorm(result, lenght);

	delete[] result;

	if (norm  < (approximation * (1 - normC)) / normC) {
		return true;
	}
	return false;
}

bool isTrueCondition(double *prevGuess, double *guess, int lenght, double approximation, double normC) {
	double *result = new double[lenght];

	CopyMatrix(result, guess, lenght, 1, 0, 0);
	SubtractMatrix(result, prevGuess, 1, lenght);

	double norm = EuclideanNorm(result, lenght);

	delete[] result;

	if (norm  < approximation / 100) {
		return true;
	}
	return false;
}

double NormC(double *A, int numOfVar, int numOfEq) {
	double sum = 0;
	for (int i = 0; i < numOfEq; i++) {
		for (int j = 0; j < numOfVar; j++) {
			double tmpValue = 0;
			if (i != j) {
				tmpValue = A[i * numOfVar + j] / A[i * numOfVar + i];
				tmpValue = powl(tmpValue, 2);
			}
			sum += tmpValue;
		}
	}
	return sqrt(sum);
}

double* Check(double *X, double *A, int numOfVariables, int numOfEqiations) {
	double *res = new double[numOfEqiations];

	for (int i = 0; i < numOfEqiations; i++) {
		res[i] = 0;
	}

	for (int i = 0; i < numOfEqiations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			res[i] += X[j] * A[i * numOfVariables + j];
		}
	}

	return res;
}

void *RelaxationMethod(double relaxationFactor, double *A, int numOfVar, int numOfEq, double *b, double *initialGuess, double approximation) {
	/* A*X = b */
	/* X = CX + g norm(C) will be required for the while loop condition */
	/* Xk = (1-w)*Xk-1 - w*(D^-1)*L*Xk + w*(D^-1)*R*Xk-1 + w*(D^-1) *b*/
	/* denote the multipliers so that they can be calculated in advance:
	w - relaxationFactor
	(1-w) denote factor1
	w*(D^-1)*L denote factor2
	w*(D^-1)*R denote factor3
	w*(D^-1)*b denote factor4
     */
	double normC = NormC(A, numOfVar, numOfEq);

	double factor1 = 1 - relaxationFactor;
	//--------------------------------------------------------------------
	double *diagA = new double[numOfVar * numOfEq];
	InitDiagMatrixFromCoefMatrix(A, diagA, numOfVar, numOfEq);
	InitInvertibleDiagMatrix(diagA, numOfVar, numOfEq);

	double *L = new double[numOfVar * numOfEq];
	InitLowerTriagMatrWithoutDiag(A, L, numOfVar, numOfEq);

	double *factor2 = MultMatrixAB(diagA, numOfVar, numOfEq, L, numOfVar, numOfEq);
	MultMatrixByNum(factor2, numOfVar, numOfEq, relaxationFactor);
	//--------------------------------------------------------------------
	double *R = new double[numOfVar * numOfEq];
	InitUpperTriagMatrWithoutDiag(A, R, numOfVar, numOfEq);

	double *factor3 = MultMatrixAB(diagA, numOfVar, numOfEq, R, numOfVar, numOfEq);
	MultMatrixByNum(factor3, numOfVar, numOfEq, relaxationFactor);
	//--------------------------------------------------------------------
	double *factor4 = MultMatrixAB(diagA, numOfVar, numOfEq, b, 1, numOfEq);
	MultMatrixByNum(factor4, 1, numOfEq, relaxationFactor);
	//--------------------------------------------------------------------

	double *prevX = new double[numOfEq];
	initMatrixByZero(prevX, 1, numOfEq);
	double tmpValue = 0;//for a separate calculation of the components of the solution vector

	do {
		CopyMatrix(prevX, initialGuess, 1, numOfEq, 0, 0);
		for (int i = 0; i < numOfEq; i++) {
			tmpValue += factor1 * initialGuess[i];

			double sum = 0;
			for (int j = 0; j < numOfVar; j++) {
				sum += factor2[i * numOfVar + j] * initialGuess[j];
			}

			//sum *= relaxationFactor;

			tmpValue -= sum;

			sum = 0;
			for (int j = 0; j < numOfVar; j++) {
				sum += factor3[i * numOfVar + j] * initialGuess[j];
			}

			//sum *= relaxationFactor;

			tmpValue -= sum;

			tmpValue += factor4[i];

			initialGuess[i] = tmpValue;
			tmpValue = 0;
		}
	} while (!isTrueCondition(prevX, initialGuess, numOfEq, approximation, normC));

	//--------------------------------------------------------------------
	delete[] diagA;
	delete[] L;
	delete[] factor2;
	delete[] R;
	delete[] factor3; 
	delete[] factor4;

	return initialGuess;
}

void PrintMatrix(double *matrix, int numOfVariables, int numOfEquations) {
	ofstream file;
	file.open("result.txt");

	file.precision(10);

	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file << matrix[i * numOfVariables + j] << "     ";
		}
		file << endl;
	}
	file << "--------------------------------------";
	file.close();
}

int main(void) {
    fileOut.open("matr3.txt");

	int numOfVariables, numOfEquations;
    GetMartixSize(&numOfVariables, &numOfEquations);

	//calculate the matrix of coefficients
	double *A = new double[numOfVariables * numOfEquations];
	GetMatrix(numOfVariables, numOfEquations, A);

	//calculate the vector of right parts
    double *b = new double[1 * numOfEquations];
	GetMatrix(1, numOfEquations, b);

	fileOut.close();

	//------------------------------------------------------

	double *X = new double[1 * numOfEquations];
	initMatrixByZero(X, 1, numOfEquations);

	double relaxationFactor = 1.5; 
	double approximation = 0.01;

	RelaxationMethod(relaxationFactor, A, numOfVariables, numOfEquations, b, X, approximation);
	//X = Check(b, A, numOfVariables, numOfEquations);
	PrintMatrix(X, 1, numOfEquations);

    //------------------------------------------------------
	delete[] A;
	delete[] b;
	delete[] X;
}