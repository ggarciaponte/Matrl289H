#include <iostream>
#include <cmath>
#include <vector>
#include "eigen-git-mirror/Eigen/Dense"
#include <cstdlib>
#include <random>
#include <unordered_map>

#define PREC 1e-6
#define kB 0.00008617
#define Vp 0.0
#define Vnn 0.030
#define V0 -0.060

///Function returns the sum of the values in a given vector-list.
double list_sum(std::vector<double>* list_ptr) 
{
    std::vector<double>& list = *list_ptr;

    double sum = accumulate(list.begin(), list.end(), 0.0);
    return sum;
}

///Function returns the average value of a given vector-list.
double average_value(std::vector<double>* list_ptr) 
{
    std::vector<double>& list = *list_ptr;

    double average = list_sum(&list) / list.size();
    return average;
}

///Function returns the square of a given value.
double value_squared(double value) 
{
    double square = pow(value, 2);
    return square;
}

///Function calculates the delta energy of a particular microstate by calculating the energy 
/// of the the most recently changed site and its nearest neighbors.
double calculate_delta_energy(Eigen::MatrixXi* lattice_ptr, int index_i, int index_j) 
{
    Eigen::MatrixXi& lattice = *lattice_ptr;

    int lattice_Nx = lattice.rows();
    int lattice_Ny = lattice.cols();
    
    ///These are the nearest neighbor indeces, considering boundary conditions.
    int xnn_right = (index_i + 1 + lattice_Nx) % lattice_Nx;
    int xnn_left = (index_i - 1 + lattice_Nx) % lattice_Nx;
    int ynn_up = (index_j + 1 + lattice_Ny) % lattice_Ny;
    int ynn_down = (index_j - 1 +lattice_Ny) % lattice_Ny;

    double delta_energy = 2 * lattice(index_i, index_j) * (Vp + Vnn * (lattice(xnn_right, index_j) + lattice(index_i, ynn_up) + lattice(xnn_left, index_j) + lattice(index_i, ynn_down)));

    return delta_energy;
}

///Function calculates the Boltzmann factor of the change in grand potential,
///exp(delta omega); delta omega (delat grand potential) = delta energy - mu * delta N. 
double calculate_boltzmann_factor(double delta_energy, double chem_pot, double delta_N, double temp) 
{
    double beta = 1 / (kB * temp);
    double delta_grand_pot = delta_energy - chem_pot * delta_N;
    double exp_term = -1 * beta * delta_grand_pot;
    double boltzmann_factor = std::exp(exp_term);
    return boltzmann_factor;
}

///Function calculates the total energy of a crystal in a particular microstate.
double calculate_total_energy(Eigen::MatrixXi* lattice_ptr)
{
    Eigen::MatrixXi& lattice = *lattice_ptr;
    
    ///Defining the dimentions of the lattice.
    int lattice_Nx = lattice.rows();
    int lattice_Ny = lattice.cols();
    
    std::vector<double> list_point_energies;
    std::vector<double> list_nearest_neighbor_energies;

    for (int i = 0; i < lattice_Nx; i++) 
    {
        for (int j = 0; j < lattice_Ny; j++) 
        {
            ///Calculates the point energy of one lattice site and pushes it back
            ///into a vector of point energyies.
            int point_energy = lattice(i, j);
            list_point_energies.push_back(point_energy);

            ///Nearest neighbor indices, considering boundary conditions.
            int right_xnn = (i + 1 + lattice_Nx) % lattice_Nx;
            int up_ynn = (j + 1 + lattice_Ny) % lattice_Ny;
            
            ///Calculates the nearest neighbor energy of one lattice site and pushes
            ///it back into a vector of nearest neighbor energies.
            int nearest_neighbor_energy = lattice(i, j) * lattice(right_xnn, j) + lattice(i, j) * lattice(i, up_ynn);
            list_nearest_neighbor_energies.push_back(nearest_neighbor_energy);
        }
    }

    ///Calculates the sum of all point energies.
    int point_energy_sum = list_sum(&list_point_energies);
    ///Calculates the sum of all nearest neighbor energies.
    int nearest_neighbor_energy_sum = list_sum(&list_nearest_neighbor_energies);

    ///Calculates the total energy.
    double total_energy = V0 + Vp * point_energy_sum + Vnn * nearest_neighbor_energy_sum;

    return total_energy;
}

///Function calculates the heat capacity from a list of grand potentials.
double calculate_heat_capacity(std::vector<double>* list_of_grand_pot_ptr, int temp) 
{
    std::vector<double>& list_of_grand_pot = *list_of_grand_pot_ptr; 
    std::vector<double> list_squared_grand_pot;
    
    ///Goes through the list of grand potentials and squares each value.
    ///Adds each squared grand potential value to a list of squared grand potentials.
    for (int i = 0; i < list_of_grand_pot.size(); i++) 
    {
        double grand_pot_squared = value_squared(list_of_grand_pot[i]);
        list_squared_grand_pot.push_back(grand_pot_squared);
    }

    ///This is the average of the squared grand potential values.
    double average_of_squared_grand_pot = average_value(&list_squared_grand_pot);
    ///Finds the average of all the (unsquared) grand potential values.
    double average_grand_pot = average_value(&list_of_grand_pot);
    ///Finds the square of the average of grand potentials.
    double average_grand_pot_squared = value_squared(average_grand_pot);
    double temp_squared = value_squared(temp);

    double heat_capacity = (average_of_squared_grand_pot - average_grand_pot_squared) / (kB * temp_squared);

    return heat_capacity;
}

bool accept_microstate(double delta_grand_pot, double boltzmann_factor, std::uniform_real_distribution<> distribution, std::mt19937_64 generator) 
{
    double rand_num_btwn_0_1 = distribution(generator);

    ///Change in grand canonical energy is less than zero.
    ///Condition 1 is met. New microstate will be accepted. 
    if (delta_grand_pot < 0) 
    {
        return true;
    }
    
    ///Change in grand canonical energy is not less than zero but has a chance of being
    ///thermally excited if the boltzmann factor is greater than a random number between 0 & 1.
    ///Condition 2 is met. New microstate will be accepted.
    else if (boltzmann_factor > rand_num_btwn_0_1) 
    {
        return true;
    }

    ///Neither condition is met. New microstate will be rejected.
    else 
    {
        return false;
    }
}

///Monte Carlo Metropolis function. Function performs an importance sampling for a given
///temperature and chemical potential.
///Function returns an average number of atoms that will occupy a crystal at a given
///temperature and chemical potential as well as a heat capacity.
std::unordered_map<std::string, double> Metropolis_MC(int temp, double chem_pot, Eigen::MatrixXi* lattice_ptr, int* N_ptr, double* total_energy_ptr) 
{
    Eigen::MatrixXi& lattice = *lattice_ptr;
    int& N = *N_ptr;
    double& total_energy = *total_energy_ptr;

    ///Generates a random decimal number in distribution 0 and 1.
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    ///Defining the dimensions of the lattice as variables.
    int lattice_Nx = lattice.rows();
    int lattice_Ny = lattice.cols();
    
    ///Number of lattice points in the lattice.
    int num_lattice_pts = lattice.size();

    ///Creating vectors to push back grand canonical energies and number of atoms in the
    ///lattice of accepted microstates.
    std::vector<double> list_accepted_grand_potentials;
    std::vector<double> list_accepted_N;
    
    int vector_size = 3000 * num_lattice_pts;

    ///Making the vectors large enough to contain all values that will be pushed back.
    list_accepted_grand_potentials.reserve(vector_size);
    list_accepted_N.reserve(vector_size);

    ///Loops through 5000 Monte Carlo passes.
    for (int m = 0; m < (5000 * num_lattice_pts); m++) 
    {
        ///Generates random indices to change a random site in the lattice.
        int rand_i = int(dis(gen) * lattice_Nx);
        int rand_j = int(dis(gen) * lattice_Ny);

        ///Changes the microstate of the crystal by either placing or taking away
        ///an atom at a randomly selected site.
        lattice(rand_i, rand_j) = lattice(rand_i, rand_j) * -1;
        
        double delta_energy = calculate_delta_energy(&lattice, rand_i, rand_j);
        int delta_N = lattice(rand_i, rand_j);

        double delta_grand_potential = delta_energy - chem_pot * delta_N;
        double boltzmann_factor = calculate_boltzmann_factor(delta_energy, chem_pot, delta_N, temp);
        
        ///For the first 1000 MC passes, the lattice is being changed, and N updated,
        ///but no energies or Ns are being stored.
        if (m < (2000 * num_lattice_pts)) 
        {
            ///If new microstate meets conditions, N is updated.
            if (accept_microstate(delta_grand_potential, boltzmann_factor, dis, gen) == true) 
            {
                N = N + delta_N;
            }

            ///If new microstate does not meet conditions, it is changed to previous microstate.
            else 
            {
                lattice(rand_i, rand_j) = lattice(rand_i, rand_j) * -1;            
            }
        }

        ///After the first 1000 passes, if new microstate meets conditions, microstate
        ///is accepted. Grand canonical energy and number of atoms get stored in vectors.
        else if (accept_microstate(delta_grand_potential, boltzmann_factor, dis, gen) == true) 
        {
            total_energy = total_energy + delta_energy;
            N = N + delta_N;
            double grand_potential = total_energy - chem_pot * N;
        
            ///Accepts new microstate.
            list_accepted_grand_potentials.push_back(grand_potential);
            list_accepted_N.push_back(N);
        }

        ///Neither condition is met. Microstate is changed to previous configuration. 
        ///Grand canonical energy and N of previous configuration added again to vectors.
        else 
        {
            lattice(rand_i, rand_j) = lattice(rand_i, rand_j) * -1;
            double grand_potential = total_energy - chem_pot * N;

            ///Accepts previous microstate once again.
            list_accepted_grand_potentials.push_back(grand_potential);
            list_accepted_N.push_back(N);
        }
    }

    double heat_capacity = calculate_heat_capacity(&list_accepted_grand_potentials, temp);
    double average_N = average_value(&list_accepted_N);

    std::unordered_map<std::string, double> sites_and_heat_capacity;

    ///Heat Capacity and N (occupies sites) are stored in a map.
    sites_and_heat_capacity["Occupied Sites"] = average_N;
    sites_and_heat_capacity["Heat Capacity"] = heat_capacity;

    return sites_and_heat_capacity;
}

int main() 
{
    ///Lattice/crystal is initialized as a matrix filled with ones.
    Eigen::MatrixXi lattice = Eigen::MatrixXi::Ones(100, 100);

    ///All ones in the lattice are converted to negative ones, so lattice is
    ///initialized as an empty lattice.
    lattice = lattice * -1;

    double num_lattice_points = lattice.size();

    ///Initializing list of temperature and chemical potential with specific ranges.
    Eigen::ArrayXd list_temp = Eigen::ArrayXd::LinSpaced(25, 100, 1000);
    Eigen::ArrayXd list_chem_pot = Eigen::ArrayXd::LinSpaced(75, -0.5, 0.5);

    ///Initializing matrices where atomic fraction and heat capacity values will be stored.
    Eigen::MatrixXd atomic_fraction_matrix(25, 75);
    Eigen::MatrixXd heat_capacity_matrix(25, 75);
    
    ///The number of atoms on the lattice is initialized as zero (since the lattice
    ///is initialized as empty).
    int N = 0;

    ///The total energy is calculated by visiting each lattice point once. Consequent
    ///total energies are calculated by adding delta energy to this initial total energy.
    double total_energy = calculate_total_energy(&lattice);
    
    ///Two for loops to loop through all chemical potentials and temperatures.
    for (int i = 0; i < list_chem_pot.size(); i++) 
    {
        for (int j = 0; j < list_temp.size(); j++) 
        {
            std::unordered_map<std::string, double> atoms_and_heat_cap = Metropolis_MC(list_temp(j), list_chem_pot(i), &lattice, &N, &total_energy);
            
            double atomic_fraction = atoms_and_heat_cap["Occupied Sites"] / num_lattice_points;
            
            ///Storing atomic fraction and heat capacity values in matrices.
            atomic_fraction_matrix(j, i) = atomic_fraction;
            heat_capacity_matrix(j, i) = atoms_and_heat_cap["Heat Capacity"];
        }
    }

    std::cout << atomic_fraction_matrix << std::endl << std::endl;
    std::cout << heat_capacity_matrix << std::endl << std::endl;

    return 0;
}         


