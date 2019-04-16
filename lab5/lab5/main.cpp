#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
const double EPS = 1e-10;

void uzupelnijTabele(double **A, double *B) {

    A[0][0] = 1.0;
    A[0][1] = 20.0;
    A[0][2] = -30.0;
    A[0][3] = -4.0;
    A[1][0] = 4.0;
    A[1][1] = 20.0;
    A[1][2] = -6.0;
    A[1][3] = 50.0;
    A[2][0] = 9.0;
    A[2][1] = -18.0;
    A[2][2] = 12.0;
    A[2][3] = -11.0;
    A[3][0] = 16.0;
    A[3][1] = -15.0;
    A[3][2] = 14.0;
    A[3][3] = 130.0;

    B[0] = 0.0;
    B[1] = 114.0;
    B[2] = -5.0;
    B[3] = 177.0;
}


void wyswietlMacierzA(double **A, int n)
{
    cout.setf(ios::fixed);
    cout.precision(2);
    cout << "Macierz A:\n";
    for (int i = 0; i<n; i++)
    {
        for (int j = 0; j<n; j++)
            cout << setw(7) << A[i][j];
        cout << endl;
    }
    cout << endl;
}

void wyswietlWektorB(double *B, int n)
{
    cout.setf(ios::fixed);
    cout.precision(2);
    for (int i = 0; i<n; i++)
        cout << setw(7) << B[i] << endl;

    cout << endl << endl;
}

void wypisz_LU(double **A, int *I, int n) {

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "Macierz U:\n";
    int i, j;
    for (i = 0; i<n; i++) {
        for (j = 0; j<n; j++) {
            if (i>j)
                cout << setw(7) << 0.00;
            else
                cout << setw(7) << A[I[i]][j];
        }
        cout << endl;
    }

    cout << endl << endl;
    cout << "Macierz L:\n";
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (j > i)
                cout << setw(7) << 0.0;
            else if(j == i)
                    cout << setw(7) << 1.0;
                else
                    cout << setw(7) << A[I[i]][j];
        }
        cout << endl;
    }
    cout << endl;
}



bool wybor_elementu_podstawowego(double **A, int *X, int i,int n) {
    int imax=1;
    for (int m = i; m < n; m++) {
        if (fabs(A[X[m]][i]) > EPS)
            imax = m;
    }

    if(imax==i)
        return false;

    // Zamianiamy wiersz macierzy

    cout << "Zamieniono elementy z wiersza " << i+1 << " z wierszem " << imax+1 << "." << endl;

    X[i] = imax;
    X[imax] = i;

    return true;
}


void dekompozycjaLU(double **A, int *X, int n) {

    int i, j, k, ii; //ii=i-1


    X[0] = 0;
    X[1] = 1;
    X[2] = 2;
    X[3] = 3;

    for (i = 1; i < n; i++)
    {
        ii = i - 1;
        if (fabs(A[X[ii]][ii]) == 0.0)
            if (wybor_elementu_podstawowego(A,X,ii,n) == false) {
                cout<<"mie mozna wybrac elementu podstawowego";
                exit(-1);
            }

        for (j = i; j < n; j++)
        {
            A[X[j]][ii]=A[X[j]][ii]/A[X[ii]][ii];
            for (k = i; k < n; ++k) {
                A[X[j]][k] = A[X[j]][k] - (A[X[ii]][k] * A[X[j]][ii]);
            }
        }

    }

}

void macierzGornaDolna(double **A, double *B, int *X, int n) {

    int i, j;
    double suma=0.0;
    //Ly=b
    //macierz trojkatna dolna:

    for (i = 0; i<n; i++){
        for (j = 0; j<i; j++)
            suma = suma + A[X[i]][j] * B[X[j]];
        B[X[i]] = B[X[i]] - suma;
        suma = 0.0;
    }


    //Ux=y
    //macierz trojkatna gorna:
    suma = 0.0;
    for (i = n-1; i >= 0; i--){
        for (j = i + 1; j<n; j++)
            suma = suma + A[X[i]][j] * B[X[j]];
        B[X[i]] = (B[X[i]] - suma) / A[X[i]][i];
        suma = 0.0;
    }


}



int main()
{

double **A, *B;
int *X, n = 4;

A = new double *[n];
for (int i = 0; i < n; i++)
    A[i] = new double[n];

B = new double[n];
X = new int[n];

uzupelnijTabele(A, B);

wyswietlMacierzA(A,n);
cout << "Wektor b:\n";
wyswietlWektorB(B, n);
dekompozycjaLU(A, X, n);
wypisz_LU(A, X, n);
macierzGornaDolna(A, B, X, n);
cout << "Wynik\n";
wyswietlWektorB(B, n);

}
