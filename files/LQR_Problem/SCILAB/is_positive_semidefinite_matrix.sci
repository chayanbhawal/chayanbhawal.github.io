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
//*********************************************************************************************************************//
// Description: This Scilab program finds if a matrix is positive semi-definite

// Input:
//   - Matrix to be checked for positive semi-definiteness
//
// Output:
//   - 1 means the matrix is positive semi-definite
//
// Assumptions to be true:
//   - The matrix is square and symmetric
//*********************************************************************************************************************//
function is_positive_semidefinite = is_positive_semidefinite_matrix(A)
    // Check if the matrix A is positive semi-definite

    // Ensure A is square
    [m, n] = size(A);
    if m ~= n
        error('Input matrix must be square.');
    end

    // Check if all eigenvalues are non-negative
    eigenvalues = (spec(A));
    no_elements = length(eigenvalues);
    is_positive_semidefinite = 1;

    for i = 1:no_elements
        if real(abs(eigenvalues(i))) >= 1D-13 && real(eigenvalues(i)) <= 0 
            is_positive_semidefinite = 0;
            break;
         end
     end
endfunction
