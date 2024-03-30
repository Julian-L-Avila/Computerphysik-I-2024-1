#include<iostream>
#include<cmath>
#include<fstream>

//Global variables
double m0=5.9736*pow(10,24);
double m1=1.9891*pow(10,30);
double e=0.016711233;
double R=149598261000;
double v=29780;
double G=6.67430*pow(10,-11);
double mu=(m0*m1)/(m0+m1);
double c=G*m0*m1;
double E=-c/(2*R);
double L=sqrt((e*e-1)*mu*c*c/(2*E));

//analitic function
double analitic(double theta){
    return ((L*L)/mu*c)/(1+cos(theta)*sqrt(1+(2*E*L*L)/(mu*c*c)));
}

int main()
{
    //Dat of analitic function
    std::ofstream datafile("analitic_funtion.dat");
    double pi=M_PI;

    for(double theta=0.0;theta<=2*3.14; theta+=0.00001)
    {
        double r_a=analitic(theta);
        datafile<<theta<<" "<<r_a<<std::endl;
    }

    datafile.close();


    //Grafic of analitic funtion
    std::ofstream scriptFile("analitic_funtion.gp");
    scriptFile<<"set term png\n";
    scriptFile<<"set output 'analitic_funtion.png'\n";
    scriptFile<<"set polar\n";
    scriptFile<<"set size square\n";
    scriptFile<<"set grid polar\n";
    scriptFile<<"x(r,theta)=r*cos(theta)\n";
    scriptFile<<"y(r,theta)=r*sin(theta)\n";
    scriptFile<<"set xlabel 'x'\n";    
    scriptFile<<"set ylabel 'y'\n";
    scriptFile<<"plot 'analitic_funtion.dat' (x($1:$2)):(y($1:$2)) with lines title 'Datos'\n";
    scriptFile.close();

    system("analitic_funtion.gp");

    return 0;


}




