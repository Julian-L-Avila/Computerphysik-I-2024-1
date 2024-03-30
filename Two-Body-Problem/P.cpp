#include<iostream>
#include<cmath>
#include<fstream>

//funcion analitica
double analitic(double theta){
    double m0;
    double m1;
    double e;
    double R;
    double v;
    double G=6.67430*pow(10,-11);
    double mu=(m0*m1)/(m0+m1);
    double L=mu*R*v;
    double c=G*m0*m1;
    double E=((e*e-1)*mu*c*c)/2*L*L;
    return ((L*L)/mu*c)/(1+cos(theta)*sqrt(1+(2*E*L*L)/(mu*c*c)));
}

int main()
{
    std::ofstream datafile("analitic_funtion.dat");
    double pi=M_PI;

    for(double theta=0.0;theta<=pi; theta+=0.00001)
    {
        datafile<<theta<<" "<<analitic(theta)<<std::endl;
        datafile<<theta<<" "<<0.0<<std::endl;
        datafile<<std::endl;
    }

    datafile.close();
}




