#include <iostream>
#include<math.h>
#include<iomanip>
#include<fstream>

#define DOKLADNOSC 1e-16
#define ITERACJE 1000

using namespace std;

fstream plik;

//funkcja sinus
double funkcja1(double x)
{
    return sin(x/4)*sin(x/4) - x;
}

//funkcja tangens
double funkcja2(double x)
{
    return tan(2*x)-x-1;
}

//liczenie pochodnej funkcji
double pochodna(double(*funkcja)(double), double x)
{
    return (funkcja(x + DOKLADNOSC) - funkcja(x)) / DOKLADNOSC;
}

double sinx_pochodna(double x) {
    return (sin(x / 2.0)) / 4.0 - 1;
}

double picard_sin(double x)
{
    return sin(x/4)*sin(x/4);
}

double picard_tan(double x)
{
    return atan(x+1)/2;
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




void metodaPicarda(double (*funkcja)(double), double (*picard_f)(double), double x0)
{
    double est=0, res=0, x;
       for(int i = 1; i<ITERACJE; i++)
    {
        if(pochodna(picard_f, x0) > 1)
        {
            cout << "Funkcja rozbiezna!" << endl;
            return;
        }
        x = picard_f(x0);
        res = fabs(funkcja(x0));
        est = fabs(x - x0);
        x0 = x;

        wypisz(i, x, est, res);
        zapisz(i, x, est, res);



        if(res < DOKLADNOSC)
        {
            plik << "residuum break" << endl;
            break;
        }
        if(est < DOKLADNOSC)
        {
            plik << "estymator break" << endl;
            break;
        }
    }
    cout << "Picard: x = " << x << endl;
}

void metodaSiecznych(double (*funkcja)(double), double a, double b)
{
    double est=0, res=0, x;
        for(int i = 1; i<ITERACJE; i++)
    {
        if(fabs(funkcja(a)) > fabs(funkcja(b)))
        {
            double temp = a;
            a = b;
            b = temp;
        }

        x = (b - a)/(funkcja(b)-funkcja(a));
        b = a;
        a = a - funkcja(a)*x;
        est = fabs(b - a);
        res = fabs(funkcja(a));

        wypisz(i, x, est, res);
        zapisz(i, x, est, res);

        if(res < DOKLADNOSC)
        {
            plik << "residuum break" << endl;
            break;
        }
        if(est < DOKLADNOSC)
        {
            plik << "estymator break" << endl;
            break;
        }
    }
    cout << "Sieczne: x = " << a << endl;
}

void metodaNewtona(double (*funkcja)(double), double x0)
{
    double est=1, res=0, x, y = 1;
    int i=0;

    while (i < ITERACJE && est > DOKLADNOSC && fabs(y) > DOKLADNOSC) {
        i++;
        x = x0;
        if(fabs(sinx_pochodna(x)) <= DOKLADNOSC){
             cout << "Rozbieżne";
             break;
        }

        x0 = x - funkcja(x) / sinx_pochodna(x);
        y = funkcja(x0);

        est = fabs(x0 - x);


         wypisz(i, x, est, res);
         zapisz(i, x, est, res);
    }
    cout << "Newton: x = " << x0 << endl;
}

void metodaBisekcji(double (*funkcja)(double), double a, double b)
{
    double est=0, res=0, x;

    if(funkcja(a)*funkcja(b)>0 || a>b)
    {
        cout << "Podano zly przedzial" << endl;
        return;
    }

    for(int i = 1; i<ITERACJE; i++)
    {
        x = (b + a)/2;
        est = fabs((b - x));
        res = fabs(funkcja(x));

        if(funkcja(a)*funkcja(x)<0)
            b = x;
        else if (funkcja(b)*funkcja(x)<0)
            a = x;
        else
        {
            cout << "Bisekcja: x = " << x << " " << (a-b)/2 << endl;
            return;
        }

        wypisz(i, x, est, res);
        zapisz(i, x, est, res);

        if(res < DOKLADNOSC) // kryterium przyblizenia wartosci funkcji blisko 0
        {
            plik << "residuum break" << endl;
            break;
        }
        if(est < DOKLADNOSC)
        {
            plik << "estymator break" << endl;
            break;
        }
    }
    cout << "Bisekcja: x = " << x << endl;
}


int main(int argc, char** argv) {
    plik.open("wyniki.txt", ios::out);

    cout << "FUNKCJA SINUS" << endl;


    metodaPicarda(funkcja1, picard_sin, -8);
    metodaBisekcji(funkcja1, -2, 4);
    metodaNewtona(funkcja1, 13);
    metodaSiecznych(funkcja1, -4, -2);


    cout << endl;
    cout << "FUNKCJA TANGENS" << endl;

    metodaPicarda(funkcja2, picard_tan, -8);
    metodaBisekcji(funkcja2, 0, 0.6);
    metodaNewtona(funkcja2, 2.5);
    metodaSiecznych(funkcja2, -5, 18);

    return 0;
}
