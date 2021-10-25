#include <iostream>
#include <armadillo>


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
    int n;
    arma::mat mat_positions;
    arma::mat mat_velocities;
    bool coulomb;

    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in, std::vector<Particle> particle_collection_in, bool coulomb_in){
        B0 = B0_in;
        V0 = V0_in;
        d = d_in;
        particle_collection = particle_collection_in;
        coulomb = coulomb_in;
        n = particle_collection.size();

        int n = particle_collection.size();
        std::cout << n << std::endl;
        mat_positions = arma::mat(3,n);
        mat_velocities = arma::mat(3,n);

        for (int i=0; i<n; i++){
          mat_positions.col(i) = particle_collection[i].position();
          mat_velocities.col(i) = particle_collection[i].velocity();
        }
    }


    // Add a particle to the trap
    void add_particle(Particle p_in){
      particle_collection.push_back(p_in);
    }

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r){
      return arma::vec({V0/(d*d)*r(0),V0/(d*d)*r(1),-2.*V0/(d*d)*r(2)});
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
      return particle_collection[i].charge()*external_E_field(mat_positions.col(i)) + particle_collection[i].charge()*arma::cross(mat_velocities.col(i),external_B_field(mat_positions.col(i)));
    }

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i, arma::mat old_pos){
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



    arma::vec acceleration(int i, arma::mat old_pos){
      return total_force(i, old_pos)/particle_collection[i].mass();
    }

    arma::mat update_velocities(arma::mat mat_velocities, arma::mat old_pos, arma::mat old_vel, double dt){
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
      for (int i=0; i<n; i++){
        particle_collection[i].position_ = mat_positions.col(i);
        particle_collection[i].velocity_ = mat_velocities.col(i);
      }
    }

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt){
      arma::mat old_pos = mat_positions;
      for (int i=0; i<n; i++){
        //Update position and velocity
        particle_collection[i].position_ = particle_collection[i].position() + dt*particle_collection[i].velocity();
        particle_collection[i].velocity_ = particle_collection[i].velocity() + dt*total_force(i,old_pos)/particle_collection[i].mass();
      }
    }

    std::vector<Particle> collection(){
      return particle_collection;
    }
};


int main(){
  //Define the time and timestep
  double time = 100;
  double dt = 0.01;
  int n = time/dt;

  //Initialize the trap with one particle
  Particle my_particle(1,40.08,arma::vec({1,0,1}),arma::vec({0,1,0}));
  Particle my_particle2(1,40.08,arma::vec({-1,0,-1}),arma::vec({0,-1,0}));

  std::vector<Particle> collection;
  collection.push_back(my_particle);
  collection.push_back(my_particle2);

  PenningTrap my_trap(96.5,9.65e8,1e4,collection,false);
  std::string filename = "test_two_part1_nocol.txt";
  std::string filename2 = "test_two_part2_nocol.txt";

  std::ofstream ofile;
  std::ofstream ffile;

  ofile.open(filename);
  ffile.open(filename2);

  ofile << collection[0].position()[0] << "     " << collection[0].position()[1] << "     " << collection[0].position()[2] << "     " << collection[0].velocity()[0] << "     "  << collection[0].velocity()[1] << "     "  << collection[0].velocity()[2] << std::endl;
  ffile << collection[1].position()[0] << "     " << collection[1].position()[1] << "     " << collection[1].position()[2] << "     " << collection[1].velocity()[0] << "     "  << collection[1].velocity()[1] << "     "  << collection[1].velocity()[2] << std::endl;

  for (int i=0; i<n; i++){
    my_trap.evolve_RK4(dt);
    std::vector<Particle> new_coll = my_trap.collection();
    ofile << new_coll[0].position()[0] << "     " << new_coll[0].position()[1] << "     " << new_coll[0].position()[2] << "     "  << new_coll[0].velocity()[0] << "     "  << new_coll[0].velocity()[1] << "     "  << new_coll[0].velocity()[2] << std::endl;
    ffile << new_coll[1].position()[0] << "     " << new_coll[1].position()[1] << "     " << new_coll[1].position()[2] << "     "  << new_coll[1].velocity()[0] << "     "  << new_coll[1].velocity()[1] << "     "  << new_coll[1].velocity()[2] << std::endl;
  }

  ofile.close();
  ffile.close();
  return 0;
}
