#include "lennardjones.h"
#include "system.h"

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

double LennardJones::sigma6()
{
    return m_sigma6;
}

double LennardJones::sigma12()
{
    return m_sigma12;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system){
    double halfLattice = system.systemSize()[0]/2.0;

    for(int i=0; i<system.atoms().size(); i++) {
       //std::cout << i <<"---" ;
       Atom *atom_i = system.atoms()[i];
       //std::cout << atom_i->mass() <<"---" ;
       for(int j=i+1; j<system.atoms().size(); j++) {
           Atom *atom_j = system.atoms()[j];
                       /* distance between atoms */
           double dx = atom_i->position[0] - atom_j->position[0];
           if (dx > halfLattice){
               dx -= system.systemSize()[0];
           }
           if (dx <= -halfLattice){
               dx += system.systemSize()[0];
           }
           double dy = atom_i->position[1] - atom_j->position[1];
           if (dy > halfLattice){
               dy -= system.systemSize()[1];
           }
           if (dy <= -halfLattice){
               dy += system.systemSize()[1];
           }
           double dz = atom_i->position[2] - atom_j->position[2];
           if (dz > halfLattice){
               dz -= system.systemSize()[2];
           }
           if (dz <= -halfLattice){
               dz += system.systemSize()[2];
           }
           double dr2 = dx*dx + dy*dy + dz*dz;
           double dr = sqrt(dr2);
           double inverse_dr2 = 1.0/dr2;
           double inverse_dr6 = inverse_dr2*inverse_dr2*inverse_dr2;
           double inverse_dr12 = inverse_dr6*inverse_dr6;
           double force = 24.0*epsilon()*(2.0*m_sigma12*inverse_dr12 - m_sigma6*inverse_dr6)*(dr/dr2);
           double fx = force*dx;
           double fy = force*dy;
           double fz = force*dz;
           atom_i->force[0] -= fx;
           atom_i->force[1] -= fy;
           atom_i->force[2] -= fz;
           atom_j->force[0] += fx;
           atom_j->force[1] += fy;
           atom_j->force[2] += fz;
           m_potentialEnergy += 4*epsilon()*((m_sigma12/inverse_dr12) - (m_sigma6/inverse_dr6));
       }
    }

    m_potentialEnergy = 0; // Remember to compute this in the loop
}
