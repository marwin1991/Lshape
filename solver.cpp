#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
using namespace std;

struct point{
    double x;
    double y;
    int neumann;
};

struct point3D{
    double x;
    double y;
    double z;
};

bool sortComperator (point3D lhs, point3D rhs)
{
    if(fabs(lhs.x - rhs.x)<1e-6) return lhs.y > rhs.y;
    return lhs.x > rhs.x;
}

int sgn(double y){
    if (y>0.0) return 1;
    if (y==0.0) return 0;
    else return -1;
}

double f(double x, double y){

    return pow(x*x,1.0/3);
    //return (pow(x*x+y*y,1.0/3)*pow(sin(acos(x/sqrt(x*x+y*y)*sgn(y))+M_PI/2),2.0/3));
    return 2*x+2*y+8;
    //return fabs(x*y);
    //return fabs( x / pow(x*x + y*y, 1.0/3) );
    //return pow((fabs(x)),2.0/3);
    //double r = sqrt(x*x+y*y);
    //double phi = acos(x/r)*sgn(y);
    //double dwie = 0.666666;
    //cout << r << " " << phi <<" "<< sin(phi+M_PI/4) <<"  "<< pow(r*sin(phi+M_PI/4),dwie)<<endl;
    //cout << (pow(r,dwie)*pow(sin(phi+M_PI/4),dwie))<<endl;
    //return (pow(r,dwie)*pow(sin(phi+M_PI/4),dwie));
}

int index (int ind, int i, int p){
    if(ind == 0) return i+2*p+1;
    if(ind == 1) return i+2*p+2;
    if(ind == 2) return i+1;
    else return i;
}
int index2 (int ind, int i, int p){
    if(ind == 0) return i+p+1;
    if(ind == 1) return i+p+2;
    if(ind == 2) return i+1;
    else return i;
}

void whitecharprint(double i){
    if(i>=0) cout << " " << i <<" | ";
    if(i<0 ) cout << i <<" | ";
}

void write_leftMatrix(double **leftMatrix, int vertax_number)
{
    cout << "    |  ";
    for ( int i = 1; i <= vertax_number; ++i )
    {
        if (i <10 )cout << " "<<i<<"   |  ";
        if (i>=10)cout << " " <<i<<"  |  ";
    }
    cout << endl << "-----";
    for ( int j = 0; j < vertax_number; ++j)
        cout << "--------";
    cout <<endl;
    for ( int i = 0; i < vertax_number; ++i )
    {
        if (i+1<10) cout << i+1 <<":  | ";
        else cout << i+1 <<": | ";
        cout << fixed;
        cout.precision(2);
        for ( int j = 0; j < vertax_number; ++j)
            whitecharprint(leftMatrix[i][j]);
        cout << endl << "-----";
        for ( int j = 0; j < vertax_number; ++j)
            cout << "--------";
        cout<<endl;
    }
}
void zeruj(point punkty[], int N){
    for (int i=0; i<N; i++)
        punkty[i].neumann=0;
}

void wypelnij_punkty1(point punkty[], int p){
    //dla dwóch pierwszych kwadratów
    double wartX = -1.0;
    double wartY = 1.0;
    double dziel = 1.0/p;
    if(fabs(dziel)<=1e-6) dziel = 0;
    int l=0;
    for (int i=0; i<(p+1); i++) {
        for (int j=0; j<(p+1)*2-1; j++) {
            punkty[l].x=wartX;
            punkty[l].y=wartY;
            wartX+=dziel;
            l++;
        }
        wartX=-1.0;
        wartY=wartY-dziel;
    }
    //trzeci kwadrat
    l=l-(p+1);
    wartX=0.0;
    wartY=0.0;
    for (int i=0; i<(p+1); i++){
        for (int j=0; j<(p+1); j++){
            punkty[l].x=wartX;
            punkty[l].y=wartY;
            wartX+=dziel;
            l++;
        }
        wartX=0.0;
        wartY=wartY-dziel;
    }
}

void wypelnij_krawedzie(point krawedzie[], point punkty[], int p){
    int l=0;
    int ind1=((p+1)*2-1)*p+1;
    int ind2=ind1-((p+1)*2-1);
    int ostatni = 3*p*p+4*p;
    //lewy bok
    for (int i=0; i<p; i++) {
        krawedzie[l].x=-1.0;
        krawedzie[l].y=(punkty[ind1].y+punkty[ind2].y)/2.0;
        punkty[ind2-1].neumann = 1;
        l++;
        ind1=ind2;
        ind2=ind2-((p+1)*2-1);
    }
    //górna krawedź
    ind1=0;
    ind2=1;
    for (int i=0; i<2*p; i++) {
        krawedzie[l].x=(punkty[ind1].x+punkty[ind2].x)/2.0;
        krawedzie[l].y=1.0;
        punkty[ind2].neumann = 1;
        l++;
        ind1++;
        ind2++;
    }
    //prawa krawedź, górna połowa
    ind2=ind1+((p+1)*2-1);
    for (int i=0; i<p; i++) {
        krawedzie[l].x=1.0;
        krawedzie[l].y=(punkty[ind1].y+punkty[ind2].y)/2.0;
        punkty[ind2].neumann = 1;
        l++;
        ind1=ind2;
        ind2=ind2+((p+1)*2-1);
    }
    //prawa krawedź, dolna połowa
    ind2=ind1+(p+1);
    for (int i=0; i<p; i++) {
        krawedzie[l].x=1.0;
        krawedzie[l].y=(punkty[ind1].y+punkty[ind2].y)/2.0;
        punkty[ind2].neumann = 1;
        l++;
        ind1=ind2;
        ind2=ind2+(p+1);
    }
    //dolna krawędź
    ind1=ostatni;
    ind2=ind1-1;
    for (int i=0; i<p; i++) {
        krawedzie[l].x=(punkty[ind1].x+punkty[ind2].x)/2.0;
        krawedzie[l].y=-1.0;
        punkty[ind1].neumann = 1;
        l++;
        ind1--;
        ind2--;
    }
}

void wypelnij_rightMatrix(double rightMatrix[], point punkty[], point krawedzie[], double phi, int N, int p){
    double dl=1.0/p*sqrt(2);
    int l=0;
    for (int i=0; i<N; i++)
        if (punkty[i].neumann==1){
            rightMatrix[i]=phi*f(krawedzie[l].x,krawedzie[l].y)*dl + phi*f(krawedzie[l+1].x,krawedzie[l+1].y)*dl;
            l++;
        }
}

void wypisz(point punkty[], int N){
    for (int i=0; i<N; i++) {
        cout << i+1 << ": " << "(" << punkty[i].x << ";" << punkty[i].y << ") " << punkty[i].neumann
        << endl;
    }
}

void wypisz_z_dl(point punkty[], int N){
    for (int i=0; i<N; i++) {
        cout << i+1 << ": " << "(" << punkty[i].x << ";" << punkty[i].y << ") " << endl;
    }
}


int main(){
    int p = 0;
    while (p<=0)
    {
        cout << "Podaj ilosc podzialu bokow kwadratu: ";
        cin >> p;
        cout << endl;
        if(p<=0)
            cout << "Niepoprawna ilosc podzialow!\n";
    }
    int squer_number = 3*p*p;
    int vertax_number = 3*p*p+ 4*p + 1;
    cout<<"Ilosc kwadratow: "<< squer_number << endl;
    cout<<"Ilosc wierzcholkow: "<< vertax_number << endl;
    //maceirz wyrazów wolnych
    double *rightMatrix = new double[vertax_number];
    point *punkty = new point[vertax_number];
    point *krawedzie = new point[p*6];
    for (int i=0; i<vertax_number; i++)
        rightMatrix[i]=0;
    zeruj(punkty, vertax_number);
    wypelnij_punkty1(punkty, p);
    wypelnij_krawedzie(krawedzie, punkty, p);
    wypelnij_rightMatrix(rightMatrix, punkty, krawedzie, 0.5, vertax_number,p);
    delete [] krawedzie;
    //macierz lewwej strony
    double **leftMatrix = new double *[vertax_number];
    for ( int i = 0; i < vertax_number; ++i )
    {
        leftMatrix[i] = new double [vertax_number];
        for ( int j = 0; j < vertax_number; ++j)
            leftMatrix[i][j]=0;
    }
    //macierz wartosci
    double M[4][4] = {2.0/3,-1.0/6,-1.0/3,-1.0/6,
        -1.0/6,2.0/3,-1.0/6,-1.0/3,
        -1.0/3,-1.0/6,2.0/3,-1.0/6,
        -1.0/6,-1.0/3,-1.0/6,2.0/3};
    //0 -> i+2p+1
    //1 -> i+2p+2
    //2 -> i+1
    //3 -> i
    //wypelnianie dla dwóch pierwszych kwadratów
    int skok = 0;
    for(int i=0; i<(2*p+1)*p - 1 ; i++)
    {
        for(int l = 0;l<4;l++)
            for(int m = 0; m<4; m++)
                leftMatrix[index(l,i,p)][index(m,i,p)]+=M[l][m];
        skok++;
        if(skok == 2*p){
            i++;
            skok=0;
        }
    }
    //wypenianie dla 3 kwadratu
    //int i=(2*p+1)*p+p; // (2*p-1)*p +p+1
    skok = 0;
    for(int i=(2*p+1)*p+p; i<(2*p+1)*(p+1)+(p+1)*p - p -2 ; i++)
    {
        for(int l = 0;l<4;l++)
            for(int m = 0; m<4; m++)
                leftMatrix[index2(l,i,p)][index2(m,i,p)]+=M[l][m];
        skok++;
        if(skok == p)
        {
            i++;
            skok=0;
        }
    }

    //zeorwanie boku dirichletai jedynka na przekatnej
    int n1 = (2*p + 1)*p; // (2p+1)*p + 1 ale pierwszy indeks to 0,0
    for (int i=0;i< p+1 ;++i)
        for ( int j = 0; j < vertax_number; ++j)
        {
            leftMatrix[n1+i][j]=0;
            if(n1+i == j) leftMatrix[n1+i][j]=1;
        }
    int n2=n1+p+1+p;
    for (int i=0;i<p; ++i )
        for ( int j = 0; j < vertax_number; ++j)
        {
            leftMatrix[n2 + (p+1)*i][j]=0;
            if(n2 + (p+1)*i == j) leftMatrix[n2 + (p+1)*i][j]=1;
        }
    //solver
    int N = vertax_number;
    // Gaussian elimination with partial pivoting
    for (int i = 0; i < N; i++) {
        // find pivot row and swap
        int max = i;
        for (int j = i + 1; j < N; j++)
            if (fabs(leftMatrix[j][i]) > fabs(leftMatrix[max][i]))
                max = j;
        //swap
        for(int k=0;k<N;k++){
        double temp;
        temp = leftMatrix[max][k];
        leftMatrix[max][k] = leftMatrix[i][k];
        leftMatrix[i][k] = temp;

        }
        //swap
        double temp2;
        temp2 = rightMatrix[max];
        rightMatrix[max] = rightMatrix[i];
        rightMatrix[i] = temp2;
        // pivot within b
        for (int j = i + 1; j < N; j++)
        {
            double m = leftMatrix[j][i] / leftMatrix[i][i];
            if(fabs(m)<1e-6) m=0;
            rightMatrix[j] -= rightMatrix[i] * leftMatrix[j][i] / leftMatrix[i][i];
        }
        // pivot within A
        for (int j = i + 1; j < N; j++) {
            double m = leftMatrix[j][i] / leftMatrix[i][i];
            if(fabs(m)<1e-6) m=0;
            for (int k = i+1; k < N; k++) {
                leftMatrix[j][k] -= leftMatrix[i][k] * m;
            }
            leftMatrix[j][i] = 0.0;
        }
    }

    double **result = new double *[vertax_number];
    for ( int i = 0; i < vertax_number; ++i )
    {
        result[i] = new double [1];
        for ( int j = 0; j < 1; ++j)
            result[i][j]=0;
    }
    for (int j = N - 1; j >= 0; j--)
    {
        double t = 0.0;
        for (int k = j + 1; k < N; k++)
            t += leftMatrix[j][k] * result[k][0];
        result[j][0] = (rightMatrix[j] - t) / leftMatrix[j][j];
    }
    for (int i = 0; i<N; i++)
        delete [] leftMatrix[i];
    delete [] leftMatrix;
    delete [] rightMatrix;
    fstream plik( "dane.dat", ios::out | ios::trunc );
    plik.close();
    point3D *tabpoint3D = new point3D [vertax_number];
    plik.open( "dane.dat", std::ios::in | std::ios::out );
    if( plik.good() == true )
    {
        for (int i =0 ; i < N ; i++)
        {
            tabpoint3D[i].x = punkty[i].x;
            tabpoint3D[i].y = punkty[i].y;
            tabpoint3D[i].z = result[i][0];
        }
        delete [] result;
        delete [] punkty;
        sort(tabpoint3D,tabpoint3D+N,sortComperator);
        double flaga = tabpoint3D[0].x;
        for (int i =0 ; i < N-1 ; i++)
        {
            if(fabs(flaga - tabpoint3D[i+1].x) > 1e-6)
            {
                flaga = tabpoint3D[i+1].x;
                plik << tabpoint3D[i].x << "  " << tabpoint3D[i].y << "  " << tabpoint3D[i].z << "\n\n";
            }
            else
            {
            plik << tabpoint3D[i].x << "  " << tabpoint3D[i].y << "  " << tabpoint3D[i].z << "\n";
            }
        }
        plik << tabpoint3D[N-1].x << "  " << tabpoint3D[N-1].y << "  " << tabpoint3D[N-1].z << "\n\n";
        plik.close();
        delete [] tabpoint3D;
    }
    system("gnuplot.exe src.gnuplot");
    return 0;
}
