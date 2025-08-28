//*********************************************************************************************************************//
// Authors: Imrul Qais (TU Delft), Chayan Bhawal (IIT Guwahati), and Debasattam Pal (IIT Bombay).
//
// Citation: Imrul Qais, Chayan Bhawal, and Debasattam Pal, "Singular LQR optimal control of DAEs: A feedback solution", Automatica, vol. 173, 2025.
// Funding Agency for the software package: Science and Engineering Research Board (SERB), Govt. of India 
// Project ID: SRG/2021/000721
// Project Duration: 31-12-2021 to 30-12-2023
//
// Created on: 29 December 2023
//*********************************************************************************************************************//
//
//**********************************************************************************************************************//
// Description: This Scilab program initializes the system matrices and cost matrices involved in the problem

// Output:
//   - LQR cost matrix represented as [Q S; S' R]
//   - System matrices E, J, L, where the DAE system is given by E(dx/dt) = Jx + Lu
// 
// Assumptions:
//   - Cost matrix should be positive semidefinite
//**********************************************************************************************************************//

function [E,J,L,LQR_cost_matrix] = input_data_matrix()
    
//  #EXAMPLE  - DAE BASED LQR PROBLEM EXAMPLE (MIMO)
  J = 1/3*[-4 2 5 8 5;2 8 14 26 20;-2 -8 -14 -23 -17;12 18 21 33 24;8 -4 -13 -19 -10];
  L = [1 1;-1 1;0 -1;1 2;-1 -2];
  E = [1 0 0 0 1;-1 -1 0 0 0;0 0 -1 -1 -1;1 2 3 3 2;-1 0 -1 -1 -2];
  Q = 1/9*[29 29 14 14 14;29 29 14 14 14;14 14 8 8 8;14 14 8 8 8;14 14 8 8 8];
  S = [0 0;0 0;0 0;0 0;0 0];
  R = [4/9 0;0 0];

//  #EXAMPLE  - DAE BASED LQR PROBLEM EXAMPLE WHERE SLOW-SPACE OF THE HAMILTONIAN IS ABSENT (SISO)
//    J = [0 0;0 1];
//    L = [1;1];
//    E = [1 0;0 0];
//    Q = [1 0;0 0];
//    S = [0;0];
//    R = 0;

//  #EXAMPLE  - DAE BASED LQR PROBLEM EXAMPLE WHERE R > 0 (SISO)
//    J = [2 1 2;1 2 1;2 1 2];
//    L = [2;2;3];
//    E = [1 -1 1;-1 0 2;0 -1 3];
//    Q = [1 1 4;1 1 4;4 4 16];
//    S = [1;1;4];
//    R = 1;
    
    
//  #EXAMPLE  - STATE-SPACE BASED SINGULAR LQR PROBLEM EXAMPLE
//    J = [3 0 -2 2 0;1 -3 2 -1 5;-2 8 3 -1 -8;-5 3 2 -2 -4;1 -5 0 0 6];
//    L = [1 0 -1 0;0 1 0 1;-1 -1 0 -2;-2 -1 1 -1;0 1 -1 2];
//    E=  eye(size(J,1),size(J,1));
//    Q = [18 -4 0 9 13;-4 15 8 -6 -5;0 8 6 -3 1;9 -6 -3 6 6;13 -5 1 6 13];
//    S = [0 0 -3 -6;0 0 9 2;0 0 3 2;0 0 -3 -4;0 0 -6 -2];
//    R = [0 0 0 0;0 0 0 0;0 0 9 0;0 0 0 4];

//  #EXAMPLE  - STATE-SPACE BASED SINGULAR LQR PROBLEM WITH THE E NON-SINGULAR BUT NOT IDENTITY
//    J = [3 0 -2 2 0;1 -3 2 -1 5;-2 8 3 -1 -8;-5 3 2 -2 -4;1 -5 0 0 6];
//    L = [1 0 -1 0;0 1 0 1;-1 -1 0 -2;-2 -1 1 -1;0 1 -1 2];
//    E=  [1 0 0 0 0;0 2 0 0 0; 0 0 3 0 0;0 0 0 1 0;0 0 0 0 1];
//    Q = [18 -4 0 9 13;-4 15 8 -6 -5;0 8 6 -3 1;9 -6 -3 6 6;13 -5 1 6 13];
//    S = [0 0 -3 -6;0 0 9 2;0 0 3 2;0 0 -3 -4;0 0 -6 -2];
//    R = [0 0 0 0;0 0 0 0;0 0 9 0;0 0 0 4];

//  #EXAMPLE  - STATE-SPACE BASED SINGULAR LQR PROBLEM
//    E = [1 0 0;0 1 0;0 0 1]; 
//    J = [1 0 1;1 0 1;1 1 0];
//    L = [0;1;0];
//    Q = [0 0 0;0 0 0;0 0 1];
//    S = [0;0;0];
//    R = 0;
//            
//  #EXAMPLE  - STATE-SPACE BASED REGULAR LQR PROBLEM EXAMPLE
//    J = [3 0 -2 2 0;1 -3 2 -1 5;-2 8 3 -1 -8;-5 3 2 -2 -4;1 -5 0 0 6];
//    L = [1 0 -1 0;0 1 0 1;-1 -1 0 -2;-2 -1 1 -1;0 1 -1 2];
//    E=  eye(size(J,1),size(J,1));
//    Q = [18 -4 0 9 13;-4 15 8 -6 -5;0 8 6 -3 1;9 -6 -3 6 6;13 -5 1 6 13];
//    S = [0 0 -3 -6;0 0 9 2;0 0 3 2;0 0 -3 -4;0 0 -6 -2];
//    R = [9 0 0 0;0 4 0 0;0 0 9 0;0 0 0 4];


    //The LQR cost matrices
    LQR_cost_matrix = [Q S;S' R];
endfunction
