//*********************************************************************************************************************//
// Authors: Imrul Qais (TU Delft), Chayan Bhawal (IIT Guwahati), and Debasattam Pal (IIT Bombay).
//
// Citation: Imrul Qais, Chayan Bhawal, and Debasattam Pal, "Singular LQR optimal control of DAEs: A feedback solution", Automatica, vol. 173, 2025.
//
// Funding Agency for the software package: Science and Engineering Research Board (SERB), Govt. of India 
// Project ID: SRG/2021/000721
// Project Duration: 31-12-2021 to 30-12-2023
//
// Created on: 29 December 2023
//*********************************************************************************************************************//
//
//**********************************************************************************************************************//
//  Description: This Scilab program finds the Markov matrix M used to find basis of fast subspace 
//
//  Input:
//  - Reduced Hamiltonian system and input matrices Ar,Br,Cr
//  - Dimension of fast subspace: nf
//
//  Output:
//  - Markov matrix M of the form [0 0 ----- CrBr;0 CrBr  ------ CrArBr;| | | | |; CrBr CrArBr -----  CrAr^(nf-d-1)Br]
//**********************************************************************************************************************//
function Markov = markov(Ar,Br,Cr,nf)
    p = [size(Cr,1)];
    d = [size(Br,2)];
    
    // Initialize an empty matrix to store intermediate results
    mat = [];
    
    // Check if nf is greater than d
    if nf > d then
        // Iterate over the past steps to construct the Markov matrix
        for iter = 0:nf-d-1
            // Append the result of Cr * Ar^iter * Br to the matrix
            mat = clean([mat;Cr*Ar^iter*Br]);
        end 
    else
        mat = [zeros(p,d)]
    end
    
    // Iterate over the remaining future steps
    for iter = 1:nf-d-1
        // Shift the matrix downward and append to the left
        mat = [shift_downward(mat,p) mat];
    end
    
    // Get the size of the resulting matrix
    [rows,cols] = size(mat);
    
    // Append zero rows at the top to complete the Markov matrix
    Markov = [mat];
endfunction

function shifted_matrix = shift_downward(mat, num_shift)
    // Get the number of rows and columns in the matrix
    [rows, cols] = size(mat);

    // Validate the number of elements to be shifted
    if num_shift >= rows
        error('Number of elements to shift should be less than the number of rows in the matrix.');
    end

    // Initialize a matrix of zeros with the same size as the input matrix
    shifted_matrix = zeros(rows, cols);

    // Shift the rows downward by the specified number and append zeros at the top
    shifted_matrix((num_shift + 1):rows, :) = mat(1:(rows - num_shift), :);
endfunction
