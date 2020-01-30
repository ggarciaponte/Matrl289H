#include <iostream>
#include <cmath>
#include <vector>
#include "eigen-git-mirror/Eigen/Dense"

#define PREC 1e-6

class Symmetry_Operation {

    public:
        Eigen::Matrix3d cart_matrix;
        Eigen::Vector3d translation;

        Symmetry_Operation(Eigen::Matrix3d input_matrix, Eigen::Vector3d input_translation) {
            Eigen::Matrix3d cart_matrix = input_matrix;
            Eigen::Vector3d translation = input_translation;
        }
};

//Checks if symmetry operation cartesian matrix, S, is unitary. S*Stranspose = Indentity.
bool is_sym_op_matrix_unitary(Eigen::Matrix3d S) {
    Eigen::Matrix3d S_product = S.transpose()*S;
    if (S_product.isIdentity(PREC) == false) {
        return false;
    }
    return true;
}

//Generates a grid within a sphere from given lattice vectors and a radius.
std::vector<Eigen::Vector3d> create_grid_points(Eigen::Matrix3d lattice, int radius){
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

//Caclulates a point group (a list of possible point group symmetry operations)
//using L primes.
std::vector<Eigen::Matrix3d> calc_point_group(Eigen::Matrix3d lattice, int radius) {
    std::vector<Eigen::Matrix3d> point_group_sym_ops;
    std::vector<Eigen::Matrix3d> L_primes = calculate_L_primes(lattice, radius);
    for (int i = 0; i < L_primes.size(); i++) {
        Eigen::Matrix3d Lp = L_primes[i];
        Eigen::Matrix3d S = Lp * lattice.inverse();
        if (is_sym_op_matrix_unitary(S) == true) {
            point_group_sym_ops.push_back(S);
        }
    }
    return point_group_sym_ops;
}

//Functor compares two matrices for equality.
struct compare_sym_op_matrices {
                        compare_sym_op_matrices(Eigen::Matrix3d Sp) : Sp(Sp) {}
                        bool operator()(Eigen::Matrix3d Sm) const {return Sp.isApprox(Sm,PREC);}

                   private:
                        Eigen::Matrix3d Sp;
};

//Checks if the list of symmetry operations is a group, by checking for closure.
bool pt_group_is_closed(std::vector<Eigen::Matrix3d> list_pt_group_sym_ops) {
    for (int s1 = 0; s1 < list_pt_group_sym_ops.size(); s1++) {
        for (int s2 = 0; s2 < list_pt_group_sym_ops.size(); s2++) {
            Eigen::Matrix3d S_prime = list_pt_group_sym_ops[s1]*list_pt_group_sym_ops[s2];
            compare_sym_op_matrices compare(S_prime);
            if (std::find_if (list_pt_group_sym_ops.begin(), list_pt_group_sym_ops.end(), compare) == list_pt_group_sym_ops.end()) {
                return false;
            }
        }
    }
    return true;
}

//Function applies a symmetry operation (matrix and translation) to a basis and returns a transformed basis. 
std::vector<Eigen::Vector3d> apply_sym_op_to_basis(Eigen::Matrix3d cart_matrix, Eigen::Vector3d translation, std::vector<Eigen::Vector3d> basis) {
    std::vector<Eigen::Vector3d> transformed_basis;
    for (int i = 0; i < basis.size(); i++) {
        Eigen::Vector3d r_prime = cart_matrix * basis[i] + translation;
        transformed_basis.push_back(r_prime);
    }
    return transformed_basis;
}

//Functor compares two vectors for equality.
struct compare_vectors {
            compare_vectors(Eigen::Vector3d V1) : V1(V1) {}
            bool operator()(Eigen::Vector3d V2) const {return V1.isApprox(V2,PREC);}

        private:
            Eigen::Vector3d V1;
};

bool compare_two_basis_sets(std::vector<Eigen::Vector3d> basis_1, std::vector<Eigen::Vector3d> basis_2) {
    if (basis_1.size() != basis_2.size()) {
        return "Error: Cannot compare two basis sets of different size.";
    }
    for (int i = 0; i < basis_1.size(); i++) {
        compare_vectors compare(basis_1[i]);
        if (std::find_if (basis_2.begin(), basis_2.end(), compare) == basis_2.end()) {
            return false;
        }
    }
    return true;
}

//Function takes in a list of point group operations and a basis for a lattice. 
//The point group operations are applied to the basis and then the transformed basis is 
//mapped onto the original basisthrough translation. Translations that map transformed basis 
//back to the original position are added to the symmetry operation (consisting of a point 
//group operation matrix and a translation), which is a factor group operation. 
//The function returns a list of factor group operations.
std::vector<Symmetry_Operation> calc_factor_group(std::vector<Eigen::Matrix3d> list_of_pt_group_ops, std::vector<Eigen::Vector3d> original_basis) {
    Eigen::Matrix3d cart_matrix;
    Eigen::Vector3d translation;
    Eigen::Matrix3d dummy_cart_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d dummy_translation = Eigen::Vector3d::Zero();
    Symmetry_Operation sym_op(cart_matrix, translation);
    sym_op.cart_matrix = dummy_cart_matrix;
    sym_op.translation = dummy_translation;
    std::vector<Symmetry_Operation> factor_group_sym_ops;
    for (int i = 0; i < list_of_pt_group_ops.size(); i++) {
        sym_op.cart_matrix = list_of_pt_group_ops[i];
        std::vector<Eigen::Vector3d> transformed_basis = apply_sym_op_to_basis(sym_op.cart_matrix, sym_op.translation, original_basis);
        if (compare_two_basis_sets(original_basis, transformed_basis) == true) {
            factor_group_sym_ops.push_back(sym_op);
        }
        for (int j = 0; j < original_basis.size(); j++) {
            for (int k = 0; k < transformed_basis.size(); k++) {
                Eigen::Vector3d translation = original_basis[j] - transformed_basis[k];
                std::vector<Eigen::Vector3d> transformed_basis_2 = apply_sym_op_to_basis(sym_op.cart_matrix, translation, original_basis);
                if (compare_two_basis_sets(original_basis, transformed_basis_2) == true) {
                    sym_op.translation = translation;
                    factor_group_sym_ops.push_back(sym_op);
                }
            }
        }
    }
    return factor_group_sym_ops;
}

int main() {

    Eigen::Matrix3d FCC_lattice;
    FCC_lattice << 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5;

    int radius = 2;
    
    std::vector<Eigen::Matrix3d> FCC_pt_group = calc_point_group(FCC_lattice, radius);

    bool FCC_pt_grp_closed = pt_group_is_closed(FCC_pt_group);
    if (FCC_pt_grp_closed != true) {
        std::cout << "Error: Point group is not closed.";
        return 0;
    }

    std::vector<Eigen::Vector3d> diamond_cubic_basis;

    Eigen::Vector3d dc_basis_r1;
    dc_basis_r1 << 0.0, 0.0, 0.0;
    Eigen::Vector3d dc_basis_r2; 
    dc_basis_r2 << 0.25, 0.25, 0.25;

    diamond_cubic_basis.push_back(dc_basis_r1);
    diamond_cubic_basis.push_back(dc_basis_r2);

    std::vector<Symmetry_Operation> diamond_cubic_factor_group = calc_factor_group(FCC_pt_group, diamond_cubic_basis);

    int factor_group_num = 0;
    std::cout << "The factor group operations for a diamond cubic structure are: \n\n";
    for (int i = 0; i < diamond_cubic_factor_group.size(); i++) {
        factor_group_num++;
        std::cout << "Factor group operation " << factor_group_num << ":\n";
        std::cout << "Point group operation: \n";
        std::cout << diamond_cubic_factor_group[i].cart_matrix << "\n\n";
        std::cout << "Translation: \n";
        std::cout << diamond_cubic_factor_group[i].translation << "\n\n\n";
    }
    
    std::cout << "There are " << factor_group_num << " operations in the diamond cubic factor group.\n\n";

    return 0;
}
