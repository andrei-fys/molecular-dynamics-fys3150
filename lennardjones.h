#ifndef LENNARDJONES_H
#define LENNARDJONES_H
#include <math.h>

class LennardJones
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;
    double m_sigma6 = pow(m_sigma, 6);
    double m_sigma12 = pow(m_sigma6,2);
    double m_potentialEnergy = 0;

public:
    LennardJones() { }
    void calculateForces(class System &system);
    double potentialEnergy() const;
    double sigma() const;
    double sigma12();
    double sigma6();
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
};
#endif
