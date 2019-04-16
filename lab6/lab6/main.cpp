#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <exception>

using namespace std;

int N = 6; //rozmiar macierzy

void wyswietlMacierz(double *U,double *D,double *L){
    double zero = 0.0;
    cout.setf(ios::fixed);
        cout.precision(2);

        cout << "Podana Macierz:\n";

        for (int i = 0; i < N; i++) {

            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    cout <<setw(7)<< D[i];
                }
                else if(j-1==i)
                {
                    cout << setw(7) << U[j-1];
                }

                else if (j+1 == i )
                {
                    cout << setw(7) << L[j];
                }

                else {
                    cout << setw(7) << zero;
                }
            }
            cout << endl;
        }
    cout << endl;
}

void wyswietlWektor(double *B){

        for (int i = 0; i<N; i++)
            cout << B[i] << endl;

    cout << endl << endl;
}

void algorytmThomasaEta(double *U,double *D,double *L,double *ETA){

    ETA[0] = D[0];
    for(int i=1; i<N;i++){
         ETA[i] = D[i] - (L[i] * U[i - 1]) / ETA[i - 1];
    }

    wyswietlWektor(ETA);
}

void algorytmThomasaR(double *L,double *ETA,double *R,double *B){
    R[0] = B[0];
    for (int i=1;i<N;i++) {
        R[i] = B[i] - (L[i] * R[i - 1] )/ETA[i - 1];
    }

    wyswietlWektor(R);
}

void rozwiazanie(double *R,double *U,double *ETA,double *X){

    X[N-1] = pow(ETA[N-1],-1) * R[N-1];
    for(int i=N-1; i>=0 ;i--){
       X[i] = (R[i] - U[i] * X[i+1]) / ETA[i];
    }

    wyswietlWektor(X);
}


int main()
{
    double *D,*U,*L,*B;
    D = new double[N];
    U = new double[N];
    L = new double[N];

    U[0] = 2.0/3.0;
    U[1] = 5.0/6.0;
    U[2] = 9.0/10.0;
    U[3] = 13.0/14.0;
    U[4] = 17.0/18.0;
    U[5] = 0.0;

    D[0] = 30.0;
    D[1] = 20.0;
    D[2] = 10.0;
    D[3] = 10.0;
    D[4] = 20.0;
    D[5] = 30.0;

    L[0] = 0.0;
    L[1] = 3.0/4.0;
    L[2] = 7.0/8.0;
    L[3] = 11.0/12.0;
    L[4] = 15.0/16.0;
    L[5] = 19.0/20;


    B = new double[N];
    B[0] = 94.0/3.0;
    B[1] = 173.0/4.0;
    B[2] = 581.0/20.0;
    B[3] = -815.0/28.0;
    B[4] = -6301.0/144.0;
    B[5] = -319.0/10.0;

    wyswietlMacierz(U,D,L);
    cout << "Wektor" << endl;
    wyswietlWektor(B);

    double *ETA;
    ETA = new double[N];
    cout << "ETA" << endl;
    algorytmThomasaEta(U,D,L,ETA);

    double *R;
    R = new double[N];
    cout << "R" << endl;
    algorytmThomasaR(L,ETA,R,B);

    double *X;
    X = new double[N];
    rozwiazanie(R,U,ETA,X);

    return 0;
}
