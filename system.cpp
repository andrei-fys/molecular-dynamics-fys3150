#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for(Atom *atom : m_atoms) {
        for (int i=0;i<3;i++){
            if (atom->position[i] > m_systemSize[i]){
                atom->position[i] -= m_systemSize[i];
            }
            if (atom->position[i] < 0.){
                atom->position[i] += m_systemSize[i];
            }
        }
    }
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum() {
    for(Atom *atom : m_atoms) {
        m_totalVelocity += atom->velocity;
    }
    m_totalVelocity /= m_atoms.size();
    for(Atom *atom : m_atoms) {
        atom->velocity -= m_totalVelocity;
    }
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for(int i=0; i<numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k<numberOfUnitCellsEachDimension; k++) {
                Atom * atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom1->position.set(i*latticeConstant,j*latticeConstant,k*latticeConstant);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);

                Atom * atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom2->position.set(i*latticeConstant+latticeConstant/2.,j*latticeConstant+latticeConstant/2.,k*latticeConstant);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom * atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom3->position.set(i*latticeConstant,j*latticeConstant+latticeConstant/2.,k*latticeConstant+latticeConstant/2.);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom * atom4 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom4->position.set(i*latticeConstant+latticeConstant/2.,j*latticeConstant,k*latticeConstant+latticeConstant/2.);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);
            }
        }
    }
    setSystemSize(vec3(numberOfUnitCellsEachDimension*latticeConstant, numberOfUnitCellsEachDimension*latticeConstant, numberOfUnitCellsEachDimension*latticeConstant));
    /*for(int i=0; i<100; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble(0, 10); // random number in the interval [0,10]
        double y = Random::nextDouble(0, 10);
        double z = Random::nextDouble(0, 10);
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10)); // Remember to set the correct system size!
    */
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);
    m_steps++;
    m_time += dt;
}
