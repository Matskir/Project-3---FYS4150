#include <iostream>
#include <armadillo>
#include <math.h>

class Particle{

  public:
    double charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

    //Constructor
    Particle(double charge_in, double mass_in, arma::vec position_in, arma::vec velocity_in){
      charge_ = charge_in;
      mass_ = mass_in;
      position_ = position_in;
      velocity_ = velocity_in;
    }

    //Method that returns the charge
    double charge(){
      return charge_;
    }

    //Method that returns the mass
    double mass(){
      return mass_;
    }

    //Method that returns the position
    arma::vec position(){
      return position_;
    }

    //Method that returns the velocity
    arma::vec velocity(){
      return velocity_;
    }
};


class PenningTrap{
  public:
    double B0;
    double V0;
    double d;
    double ke = 1.38935333e5;
    std::vector<Particle> particle_collection;
    arma::mat mat_positions;
    arma::mat mat_velocities;
    bool coulomb;
    double f;
    double wf;
    double t;


    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_collection_in, bool coulomb_in, double f_in, double wf_in){
        B0 = B0_in;
        V0 = V0_in;
        d = d_in;
        particle_collection = particle_collection_in;
        coulomb = coulomb_in;
        f = f_in;
        wf = wf_in;
        t = 0;
    }


    // Add a particle to the trap
    void add_particle(Particle p_in){
      particle_collection.push_back(p_in);
    }

    //Method that initializes the collection using the newly added particles
    void Initialize_collection(){
      int n = particle_collection.size();
      mat_positions = arma::mat(3,n);
      mat_velocities = arma::mat(3,n);

      for (int i=0; i<n; i++){
        mat_positions.col(i) = particle_collection[i].position();
        mat_velocities.col(i) = particle_collection[i].velocity();
      }
    }

    // External time-dependent electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r){
      double V0_new = V0*(1 + f*cos(wf*t));
      return arma::vec({V0_new/(d*d)*r(0),V0_new/(d*d)*r(1),-2.*V0_new/(d*d)*r(2)});
    }

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r){
      return arma::vec({0,0,B0});
    }

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j, arma::mat old_pos){
      return ke*particle_collection[i].charge()*particle_collection[j].charge()*(mat_positions.col(i)-old_pos.col(j))/(std::pow(arma::norm(mat_positions.col(i)-old_pos.col(j)),3));
    }

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i){
      arma::vec force(3,arma::fill::zeros);
      if (arma::norm(mat_positions.col(i))>d){
        return force;
      }
      else{
        return particle_collection[i].charge()*external_E_field(mat_positions.col(i)) + particle_collection[i].charge()*arma::cross(mat_velocities.col(i),external_B_field(mat_positions.col(i)));
      }
    }

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i, arma::mat old_pos){
      int n = particle_collection.size();
      arma::vec force(3,arma::fill::zeros);
      if (coulomb){
        for (int j=0; j<n; j++){
          if (i!=j){
            force = force + force_particle(i,j, old_pos);
          }
        }
      }
      return force;
    }

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i, arma::mat old_pos){
      return total_force_external(i) + total_force_particles(i, old_pos);
    }


    //Method that computes a particles' acceleration
    arma::vec acceleration(int i, arma::mat old_pos){
      return total_force(i, old_pos)/particle_collection[i].mass();
    }

    //Method that updates the velocities of the particles to be used in RK4
    arma::mat update_velocities(arma::mat mat_velocities, arma::mat old_pos, arma::mat old_vel, double dt){
      int n = particle_collection.size();
      for (int i=0; i<n; i++){
        mat_velocities.col(i) = old_vel.col(i) + dt/2.*acceleration(i, old_pos);
      }
      return mat_velocities;
    }

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt){
      arma::mat old_pos = mat_positions;
      arma::mat old_vel = mat_velocities;

      //First guess
      arma::mat k1 = dt*mat_velocities;
      mat_positions = old_pos + k1/2.;
      mat_velocities = update_velocities(mat_velocities, old_pos, old_vel, dt);

      //Second guess
      arma::mat k2 = dt*mat_velocities;
      mat_positions = old_pos + k2/2.;
      mat_velocities = update_velocities(mat_velocities, old_pos, old_vel, dt);

      //Third guess
      arma::mat k3 = dt*mat_velocities;
      mat_positions = old_pos + k3;
      mat_velocities = update_velocities(mat_velocities, old_pos, old_vel, dt);

      //Fourth and final guess
      arma::mat k4 = dt*mat_velocities;
      mat_positions = old_pos + 1.0/6.0*(k1+2.*k2+2.*k3+k4);
      mat_velocities = update_velocities(mat_velocities, old_pos, old_vel, 2.*dt);

      //Update the positions and velocities in particle_collection
      int n = particle_collection.size();
      for (int i=0; i<n; i++){
        int n = particle_collection.size();
        particle_collection[i].position_ = mat_positions.col(i);
        particle_collection[i].velocity_ = mat_velocities.col(i);
      }
      t = t + dt;
    }

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt){
      int n = particle_collection.size();
      arma::mat old_pos = mat_positions;
      for (int i=0; i<n; i++){
        //Update position and velocity
        particle_collection[i].position_ = particle_collection[i].position() + dt*particle_collection[i].velocity();
        particle_collection[i].velocity_ = particle_collection[i].velocity() + dt*total_force(i,old_pos)/particle_collection[i].mass();
      }
    }

    //Method that returns the collection of particles. Used for printing to file
    std::vector<Particle> collection(){
      return particle_collection;
    }

    //Method that returns the number of particles left inside the Penning Trap
    int particles_inside(){
      int n = particle_collection.size();
      int numb = 0;
      for (int i=0; i<n; i++){
        if (arma::norm(particle_collection[i].position())<d){
          numb = numb + 1;
        }
      }
      //std::cout << numb << std::endl;
      return numb;

    }
};


int main(){


  //Define the time and timestep
  double time = 500;
  double dt = 0.01;
  int n = time/dt;

  //Number of particles
  int numb = 100;

  std::vector<double> numb_vector;
  std::vector<double> wf_vector;

  //While loop that runs over the frequencies
  double wf = 0.2;
  while (wf <= 0.6){

    //Create an empty PenningTrap
    std::vector<Particle> collection;
    PenningTrap my_trap(96.5,241213.1395,500,collection,false,0.4,wf);

    //Using a set seed
    arma::arma_rng::set_seed(1999);

    //A for loop that fills the collection with 100 random particles
    for (int i=0; i<numb; i++){
      arma::vec r = arma::vec(3).randn()*0.1*my_trap.d;
      arma::vec v = arma::vec(3).randn()*0.1*my_trap.d;
      Particle p(1,40.08,r,v);
      my_trap.add_particle(p);
    }

    //Initialize the collection
    my_trap.Initialize_collection();

    //Evolve in time using RK4
    for (int i=0; i<n; i++){
      my_trap.evolve_RK4(dt);
    }

    //Add the number of remaining particles to an array
    double number = my_trap.particles_inside();
    numb_vector.push_back(number);
    wf_vector.push_back(wf);

    //Increase the frequenc7
    wf += 0.02;
  }

  //Write to file
  int length = numb_vector.size();
  std::ofstream myfile;
  myfile.open ("nocoul_final_test.txt");
  for (int i=0; i<length; i++){
    myfile << numb_vector[i] << "    " << wf_vector[i] << std::endl;
  }
  myfile.close();

  return 0;
}
