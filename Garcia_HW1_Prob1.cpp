#include <iostream>
#include <cmath>
#include <vector>
#include "eigen-git-mirror/Eigen/Dense"

#define PREC 1e-6

//Checks if point group operation matrix, S, is unitary. S*Stranspose = Identity.
bool is_point_group_op_unitary(Eigen::Matrix3d S) {
    Eigen::Matrix3d S_product = S.transpose()*S;
    if (S_product.isIdentity(PREC) == false) {
        return false;
    }
    return true;
}

//Generates a grid within a sphere from given lattice vectors and a radius.
std::vector<Eigen::Vector3d> create_grid_points(Eigen::Matrix3d lattice, int radius) {
    std::vector<Eigen::Vector3d> list_grid_points;
    for (int n = -radius; n < (radius+1); n++) {
        for (int m = -radius; m < (radius+1); m++) {
           for (int k = -radius; k < (radius+1); k++) {
             Eigen::Vector3d point = n * lattice.col(0) + m * lattice.col(1) + k * lattice.col(2);
             list_grid_points.push_back(point);
           }
        }
    }
    return list_grid_points;
}

//Calculates a list of L primes (possible transformed lattiice vectors) 
//based on the newly created lattice grid.
std::vector<Eigen::Matrix3d> calculate_L_primes(Eigen::Matrix3d lattice, int radius) {
    std::vector<Eigen::Matrix3d> list_L_primes;
    std::vector<Eigen::Vector3d> points = create_grid_points(lattice, radius);
    for (int p1 = 0; p1 < points.size(); p1++) {
        for (int p2 = 0; p2 < points.size(); p2++) {
            for (int p3 = 0; p3 < points.size(); p3++) {
                Eigen::Matrix3d Lprime;
                Lprime.col(0) << points[p1];
                Lprime.col(1) << points[p2];
                Lprime.col(2) << points[p3];
                list_L_primes.push_back(Lprime);
            }
        }
    }
    return list_L_primes;
}

//Caclulates a point group (a list of possible symmetry operations)
//using L primes.
std::vector<Eigen::Matrix3d> calc_point_group(Eigen::Matrix3d lattice, int radius) {
    std::vector<Eigen::Matrix3d> list_point_group_ops;
    std::vector<Eigen::Matrix3d> L_primes = calculate_L_primes(lattice, radius);
    for (int i = 0; i < L_primes.size(); i++) {
        Eigen::Matrix3d Lp = L_primes[i];
        Eigen::Matrix3d S = Lp * lattice.inverse();
        if (is_point_group_op_unitary(S) == true) {
            
            list_point_group_ops.push_back(S);
        }
    }
    return list_point_group_ops;
}

//Functor compares two matrices for equality.
struct compare_sym_op_matrices {
                       compare_sym_op_matrices(Eigen::Matrix3d Sp) : Sp(Sp) {}
                       bool operator()(Eigen::Matrix3d Sm) const {return Sp.isApprox(Sm,PREC);}

                   private:
                       Eigen::Matrix3d Sp;
};

//Checks if the list of symmetry operations is a group, by checking for closure.
bool group_is_closed(std::vector<Eigen::Matrix3d> list_point_group_ops) {
    for (int s1 = 0; s1 < list_point_group_ops.size(); s1++) {
        for (int s2 = 0; s2 < list_point_group_ops.size(); s2++) {
            Eigen::Matrix3d S_prime = list_point_group_ops[s1]*list_point_group_ops[s2];
            compare_sym_op_matrices compare(S_prime);
            if (std::find_if (list_point_group_ops.begin(), list_point_group_ops.end(), compare) == list_point_group_ops.end()) {
                return false;
            }
        }
    }
    return true;
}

//Checks for the difference between two values.
bool almost_equal(double LHS, double RHS) {
    if (std::abs(LHS - RHS) < PREC) {
        return true;
    }
    else {
        return false;
    }
}

//Calculates eigenvalues and eigenvectors of the symp_op. Returns list of real eigenvectors of
//eigenvalues equal to 1.
std::vector<Eigen::Vector3d> eigenvectors_with_unit_eigenvals(Eigen::Matrix3d point_group_op) {
    Eigen::EigenSolver<Eigen::Matrix3d> solver(point_group_op, true);
    Eigen::Vector3<std::complex<double>> eigenvals = solver.eigenvalues();
    std::vector<Eigen::Vector3d> output_eigenvectors;
    Eigen::Matrix3<std::complex<double>> eigenvectors =solver.eigenvectors();

    int count = 0;
    for (int i = 0; i < 3; i++) {
        auto eigenval = eigenvals(i);
        if (almost_equal(eigenval.real(), 1) == true) {
            //std::cout << eigenvectors.col(i) << std::endl;
            Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
            output_eigenvectors.push_back(real_eigenvector);
            count++;
        }
    }
    if (count == 0) {
        for (int i = 0; i < 3; i++) {
            auto eigenval = eigenvals(i);
            if (almost_equal(eigenval.real(), -1) == true) {
                //std::cout << eigenvectors.col(i) << std::endl;
                Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
                output_eigenvectors.push_back(real_eigenvector);
                count++;
            }
        }
    }
    if(count == 0) {
            std::cout << "There is an error with the eigenval function.\n\n";
    }
    return output_eigenvectors;
}

//Takes in a sym op and categorizes the type of operation based on matrix properties.
std::string categorize_point_group_op(Eigen::Matrix3d point_group_op, Eigen::Matrix3d lattice) {
    double trace = point_group_op.trace();
    std::string type;
    double det = point_group_op.determinant();

    if (3 - trace < PREC) {
        type = "Identity";
        return type;
    }
    if (trace + 3 < PREC) {
        type = "Inversion";
        return type;
    }
    std::vector<Eigen::Vector3d> eigenvectors = eigenvectors_with_unit_eigenvals(point_group_op);
    if (eigenvectors.size() == 2) {
        type = "Mirror";
        return type;
    }
    else if (eigenvectors.size() == 1) {
        if (abs(det - 1) < PREC) {
            type = "Rotation";
            return type;
        }
        else if (abs(det + 1) < PREC) {
            type = "Improper Rotation";
            return type;
        }
    }
    else {
        type = "Error: Type not identified.";
        return type;
    }
}

int main() {

    Eigen::Matrix3d FCC_lattice;
    FCC_lattice << 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5;

    int radius = 2;
    
    std::vector<Eigen::Matrix3d> FCC_pt_group = calc_point_group(FCC_lattice, radius);

    bool FCC_pt_grp_closed = group_is_closed(FCC_pt_group);
    if (FCC_pt_grp_closed == true) {
        std::cout << "The point group is closed.\n\n";
    }
    else {
        std::cout << "Error: The point group is not closed.\n\n";
        return 0;
    }

    int sym_op_num = 0;
    std::cout << "The symmetry operations are: \n\n";
    for (int i=0; i<FCC_pt_group.size(); i++) {
        sym_op_num++;
        std::string sym_op_type = categorize_point_group_op(FCC_pt_group[i], FCC_lattice);
        std::cout <<"Point group operation " << sym_op_num << " is: " << sym_op_type << ".\n"; 
        std::cout << FCC_pt_group[i] << "\n\n";
    }

    
    std::cout << "There are " << sym_op_num << " operations in the FCC lattice point group.\n\n";
    

    return 0;
}
