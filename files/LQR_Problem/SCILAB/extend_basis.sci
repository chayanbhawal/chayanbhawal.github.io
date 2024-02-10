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
//  Description: This Scilab program finds extends a basis of a subspace B to a basis of another subspace W where B is contained within W 
//
//  Input:
//  - Matrix of basis of subspace B - B_basis
//  - Matrix of subspace of basis W - W_basis
//
//  Output:
//  - Extended basis such that img[B_basis new_columns] = img(W_basis)
//**********************************************************************************************************************//

function extended_basis = extend_basis(B_basis, W_basis)
    n = size(W_basis, 2); // Dimension of the vector space R^n
    d = size(B_basis, 2); // Dimension of the initial basis B

    C = B_basis; // Initialize C with the current basis B
    
    while size(C, 2) < n
        // Choose a vector v from W_basis that is not in the span of C
        v = find_new_vector(C, W_basis);
        
        // Add v to C
        C = [C, v];
    end

    extended_basis = C; // Return the extended basis for W
endfunction

function new_vector = find_new_vector(C, W_basis)
    // Find a vector in W_basis that is not in the span of C
    for i = 1:size(W_basis, 2)
        vector_candidate = W_basis(:, i);
        if ~is_in_span(vector_candidate, C)
            new_vector = vector_candidate;
            return;
        end
    end
endfunction

function result = is_in_span(v, basis)
    // Check if vector v is in the span of the basis
    result = rank([basis, v]) == rank(basis);
endfunction

