#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <cmath>
using namespace std;

typedef double(*funkcja)(double);
fstream plik;

const int N_MAX = 1000;
const double TOLX = 1e-16; //|Xn-X(n-1)|=Epsilon <=tolx(~~0)
const double TOLF = 1e-16; //|f(x)|<=tolf (~~0)
//(miejsce zerowe:
//tan w x~~2.20497296033312)
//sin x=0

//przedzia³y bisekcji:
//: sin:
const double a0_sin = -9;
const double b0_sin = 10;
const double x_0 = 2.0;
const double x_1 = 0.4;

//tan:
const double a0_tg = 2.18;
const double b0_tg = 2.23;

const double x0_tan = 2.33;
const double x1_tan = 2.25;

double sinx(double x);
double sinx_pochodna(double x);
double tanx(double x);
double tanx_pochodna(double x);
double picard_sinx(double x);
double picard_sinx_poch(double x);

double picard_TANx_POCH(double x);
void zbieznosc(funkcja f);
void wypisz(int i, double x, double y, double estymator);
void zapisz(int i, double x, double y, double estymator);

void Metoda_siecznych(funkcja f, funkcja fp, double x0, double x1);
void Metoda_Newtona(funkcja f, funkcja fp, double x1);
void Bisekcja(funkcja f, double an, double bn);
void Picard_sin(funkcja f, double x);

//***********************************************************************************************************
int main()
{

    plik.open("plik.txt", ios::out);

    cout <<"SINUS" <<endl;
    zbieznosc(picard_sinx_poch);
    Picard_sin(picard_sinx,x_1);
    Bisekcja(sinx, a0_sin, b0_sin);
    cout << endl;
    Metoda_Newtona(sinx, sinx_pochodna, x_1);
    cout << endl;
    Metoda_siecznych(sinx, sinx_pochodna, x_0, x_1);
    cout << endl;

    cout << "TANGENS" << endl;

    zbieznosc(picard_TANx_POCH);
    Bisekcja(tanx,a0_tg,b0_tg);
    cout << endl;
    Metoda_Newtona(tanx, tanx_pochodna, x1_tan);
    cout << endl;
    Metoda_siecznych(tanx, tanx_pochodna, x0_tan, x1_tan);
    cout << endl;


    //cout<<endl;

    plik.close();

    system("pause");
    return 0;
}




//********************************************************************************
void Metoda_siecznych(funkcja f, funkcja fp, double x0, double x1) {
    double x, y = 1.0, yp, estymator = 1.0;
    int    i = 0;

    cout << "Metoda siecznych\nn:     x:               residuum:        estymator: " << endl;
    plik << "Metoda siecznych\nn:     x:               residuum:        estymator: " << endl;

    while (i < N_MAX && estymator > TOLX && fabs(y) > TOLF) {
        //x=xn+2 , x1=xn+1,  x0=xn
        i++;
        //y = f(x1);
        x = x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0));
        y = f(x);

        estymator = fabs(x - x1);

        x0 = x1;
        x1 = x;


        wypisz(i, x, y, estymator);
        zapisz(i, x, y, estymator);


    }
    cout << endl<<"\nend\n\n";

}

//*************************************************************
void Metoda_Newtona(funkcja f, funkcja fp, double x1) {
    //x=Xn fp=pochodna
    //x1=Xn+1
    double x, y = 1.0, estymator = 1.0;
    int    i = 0;

    cout << "Metoda Newtona \nn:     x:               residuum:        estymator: " << endl;
    plik << "Metoda Newtona \nn:     x:               residuum:        estymator: " << endl;
    while (i < N_MAX && estymator > TOLX && fabs(y) > TOLF) {
        i++;

        x = x1;
        if (fabs(fp(x)) <= TOLF) {
            cout << "rozbie¿ne"; break;
        }

        x1 = x - f(x) / fp(x);
        y = f(x1);



        estymator = fabs(x1 - x);



        wypisz(i, x1, y, estymator);
        zapisz(i, x1, y, estymator);


    }
    cout << endl;
}


//******************************************************************
void Bisekcja(funkcja f, double an, double bn) {
    int i = 0;
    double xs = 1, y = 1, estymator, a, b;

    cout << "Metoda bisekcji \nn:     x:               residuum:        estymator: " << endl;
    plik << "Metoda bisekcji \nn:     x:               residuum:        estymator: " << endl;

    estymator = (bn - an) / 2.0;


    while (i < N_MAX && estymator>TOLX && fabs(y)>TOLF)
    {
        xs = (bn + an) / 2.0;//xs-srodek przedzia³u
        i++;
        y = f(xs);
        a = f(an);
        b = f(bn);


        if (a*b > 0) {
            cout << "Przedzial nie spelnia zalozen" << endl;
            break;
        }

        if (y*a < 0)
            bn = xs;//[a:s]

        if (y*b < 0)
            an = xs;//[b:s]

        wypisz(i, xs, y, estymator);
        zapisz(i, xs, y, estymator);


        estymator = (bn - an) / 2.0;
    }
    cout << endl;
}

//*******************************************************************

void Picard_sin(funkcja f, double x) {
    int i = 0;
    double y = 1, estymator = 1;


    cout << "Picard \nn:     x:               residuum:        estymator: " << endl;
    plik << "Picard \nn:     x:               residuum:        estymator: " << endl;



    while (i < N_MAX && estymator>TOLX && fabs(y)>TOLF) {
        i++;
        y = picard_sinx(x);
        estymator = fabs(x - y);
        x = y;//prosta y=x

        wypisz(i, x, y, estymator);
        zapisz(i, x, y, estymator);

    }

    cout << endl;
}


//***********************************************************************************************************
double sinx(double x) {
    return sin(x / 4.0)*sin(x / 4.0) - x;
}

double sinx_pochodna(double x) {
    return (sin(x / 2.0)) / 4.0 - 1;

}

double tanx(double x) {
    return tan(2 * x) - x - 1;
}

double tanx_pochodna(double x) {
    return 2 / (cos(2 * x)*cos(2 * x)) - 1;
}

double picard_sinx(double x) {
    return sin(x / 4.0)*sin(x / 4.0);
}


double picard_TANx_POCH( double x) {
    return 2.0 / (cos(2 * x)*cos(2 * x));
}

double picard_sinx_poch(double x) {
    return (1.0 / 4.0)*sin(x / 2);
}



void wypisz(int i,double x,double y, double estymator) {

    cout.width(3);
    cout << internal << i << "  ";
    cout.precision(7);
    cout.width(15);
    cout << internal << scientific << x << "  ";
    cout.width(15);
    cout << internal << fabs(y) << "  ";
    cout.width(15);
    cout << internal << estymator << "  " << endl;
}

void zapisz(int i, double x, double y, double estymator)
{
    if (plik.good()) {
        plik.width(3);
        plik << left << right << i << "  ";
        plik.precision(7);
        plik.width(15);
        plik << internal << scientific << x << "  ";
        plik.width(15);
        plik << internal << fabs(y) << "  ";
        plik.width(15);
        plik << internal << estymator << "  " << endl;
    }

}



void zbieznosc(funkcja f) {
    double x = -3.14;
    for (x;x<3.14;x=x+0.1){

    if (f(x) > 1) {
        cout << "Picard: metoda nie jest zbiezna"<<endl<<endl;
        break;
    }

    }
}







