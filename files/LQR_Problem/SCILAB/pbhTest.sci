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
// Description: This Scilab program is used to check if an LTI system (A,B) is stabilizable
//
// Input:
//   - System matrices A, B, where the state-space system is given by (dx/dt) = Ax + Bu
//
// Output:
//   - isStabilizable = 0, if system is stabilizable. Else 1.
//**********************************************************************************************************************//

function isStabilizable = pbhTest(A, B)
    // Perform Popov-Belevitch-Hautus (PBH) test for stabilizability
    
    // Number of states
    n = size(A, 1);
    
    // Eigenvalues of A
    eigenvalues_A = spec(A);
    
    isStabilizable = 0;//True
    
    for i = 1:length(eigenvalues_A)
        s = eigenvalues_A(i);
        
        // Check the rank condition for each eigenvalue
        if rank([s*eye(n,n) - A, B]) ~= n
            isStabilizable = 1;//False
            break;
        end
    end
endfunction
