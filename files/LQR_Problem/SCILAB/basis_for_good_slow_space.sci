//*********************************************************************************************************************//
// Authors: Imrul Qais (TU Delft), Chayan Bhawal (IIT Guwahati), and Debasattam Pal (IIT Bombay).
//
// Citation: Imrul Qais, Chayan Bhawal, and Debasattam Pal, "Singular LQR optimal control of DAEs: A feedback solution", Under Review, 2023.
//
// Funding Agency for the software package: Science and Engineering Research Board (SERB), Govt. of India 
// Project ID: SRG/2021/000721
// Project Duration: 31-12-2021 to 30-12-2023
//
// Created on: 29 December 2023
//*********************************************************************************************************************//
//
//**********************************************************************************************************************//
//  Description: This Scilab program finds basis vectors for good slow subspace of an output-nulling Hamiltonian system. 
//
//  Input:
//  - Reduced Hamiltonian matrix pair: (Er,Hr)  
//  - Dimension of good slow subspace: ns
//
//  Output:
//  - A matrix Vlambda whose columns are a basis of  good slow subspace of a output-nulling Hamiltonian system 
//**********************************************************************************************************************//
function [Vlambda,Lambda_set] = basis_for_good_slow_space(Er,Hr,ns)
    // Compute the eigenvalues of reduced Hamiltonian matrix pair (Er,Hr)
    eig = roots(det(s*Er-Hr));
    
    //Construction of Lambda-set
    
    //Passing the vector of eigen-values of (Er,Hr) to the function to get the vector of eigenvalues in -ve complex plane
    Lambda_set = find_negative_elements(eig);    
    
    // Initialize a matrix to store basis vectors for the slow subspace
    Vlambda = [];
    Gamma = []
    // Iterate over the selected number of Lambda-set elements to find a basis of the corresponding eigenspace 
    for iter = 1:ns
        // Append the kernel (null space) of (Lambda_set(iter)*Er - Hr) to Vlambda
        Vlambda = [Vlambda, kernel(Lambda_set(iter)*Er-Hr)];
        Gamma = blockdiag(Gamma,Lambda_set(iter))
    end
 endfunction 
 
//Find elements in a vector that have negative real part

 function negative_elements = find_negative_elements(vector)
    // Initialize an empty array to store negative elements
    negative_elements = [];

    // Loop through each element in the vector
    for i = 1:length(vector)
        // Check if the current element is negative
        if real(vector(i)) < 0
            // If negative, add it to the array
            negative_elements = [negative_elements, vector(i)];
        end
    end
endfunction
 
