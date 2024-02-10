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
// Description: This Scilab program finds feedback gain matrices to solve the LQR problem for DAE LTI systems irrespective of the index of the system

// Input:
//   - LQR cost matrix represented as [Q S; S' R] 
//   - System matrices E, J, L, where the state-space LTI system is given by E(dx/dt) = Jx + Lu
//   - No. of states 
//   - No. of inputs 
//
// Output:
//   - The state-feedback gain matrices Fp and Fd such that the control law u = Fp x + Fd (dx/dt) solves the LQR problem
//   - The closed loop system matrices of the system after u = Fp x + Fd (dx/dt) is employed. 
//   - The closed loop system is Ecl (dx/dt) = Acl x 
//
// Assumptions to be true:
//   - Cost matrix should be positive semidefinite where the input cost matrix R may or may not be singular
//   - The system must be impulse controllable
//*********************************************************************************************************************//
function [Fp,Fd,ns,nf,expect_cl_poles] = LQR_problem_DAE(E,J,L,LQR_cost_matrix,no_states,no_inputs)  
    //Dimension of the slow subsystem of the DAE
    dim_slow_system = rank(E);
    //Dimension of the fast subsystem of the DAE
    dim_fast_system = no_states - dim_slow_system;
    
    //Computation of matrices to make E block diagonal Phat*E*That = [I 0;0 0]
    [Phat,That] = block_Ematrix(E);
    
    //Changing the input space basis uhat = That X u and post multiplying the system with Phat
    
    //The system matrix J in the new input basis
    Jhat = Phat*J*That;
    //Partitioning the system matrix conforming to the slow and fast subsystem partition
    J11 = Jhat(1:dim_slow_system,1:dim_slow_system);
    J12 = Jhat(1:dim_slow_system,dim_slow_system+1:no_states);
    J21 = Jhat(dim_slow_system+1:no_states,1:dim_slow_system);
    J22 = Jhat(dim_slow_system+1:no_states,dim_slow_system+1:no_states);
    
    //The input matrix L in the new input basis
    Lhat = Phat*L;
    //Partitioning the system matrix conforming to the slow and fast subsystem partition
    B1 = Lhat(1:dim_slow_system,:);
    B2 = Lhat(dim_slow_system+1:no_states,:);
    
    //Check if the system is impulse controllable and if yes continue with solving the problem, else return empty matrices
    if rank([J22 B2]) == dim_fast_system then
        disp('System is impulse controllable');
        
        //Choice of F2 - if J22 is full rank, then F2 = 0, else a random F2 is chosen 
        if rank(J22) == dim_fast_system then
            F2 = zeros(no_inputs,dim_fast_system);
        else
            F2 = rand(no_inputs,dim_fast_system);
        end
        
        //Converting the DAE based LQR problem to a state-space based LQR problem using a feedback control law
        
        //Feedback matrix to convert the system to an intermediate form using the control law uhat = F xhat + u
        F = [zeros(no_inputs,dim_slow_system) F2]*inv(That);
        
        //Construction of matrices to convert the intermediate system to Weistrass canonical form
        
        //PET = [I 0;0 0] and P(J+LF)T = [A 0;0 I]
        P = [eye(dim_slow_system,dim_slow_system) -(J12 + B1*F2)*inv(J22 + B2*F2);
            zeros(dim_fast_system,dim_slow_system) eye(dim_fast_system,dim_fast_system)]*Phat;
        T = That*[eye(dim_slow_system,dim_slow_system) zeros(dim_slow_system,dim_fast_system);
                  -inv(J22 + B2*F2)*J21 inv(J22+B2*F2)];
        
        //The slow subsytem matrices
        A = (P*(J+L*F)*T)(1:dim_slow_system,1:dim_slow_system);
        B = (P*L)(1:dim_slow_system,:);

        //Construction of transformation matrix to change the cost matrix to the new state space and new input basis
        Tran = [T zeros(no_states,no_inputs);
               F*T eye(no_inputs,no_inputs)]*[eye(dim_slow_system,dim_slow_system) zeros(dim_slow_system,no_inputs);                     zeros(dim_fast_system,dim_slow_system) -B2;zeros(no_inputs,dim_slow_system) eye(no_inputs,no_inputs)];
        
        //Transformed cost matrix 
        Cost_trf = Tran'*LQR_cost_matrix*Tran;
        
        //Partitioning the cost matrix conforming to state-space and input-space dimension
        Qhat = clean(Cost_trf(1:dim_slow_system,1:dim_slow_system));
        Shat = clean(Cost_trf(1:dim_slow_system,dim_slow_system+1:dim_slow_system+no_inputs));
        Rhat = clean(Cost_trf(dim_slow_system+1:dim_slow_system+no_inputs,dim_slow_system+1:dim_slow_system+no_inputs));
     else
         //System is not impulse controllable and hence cannot be solved for arbritrary initial condition
         disp('Problem not solvable for arbitrary initial condition - System not impulse controllable');
         A = [];
         B = [];
         Qhat = [];
         Shat = [];
         Rhat = []; 
     end
     
     //Executing function to solve LQR problem for state-space (A,B) and cost matrix [Qhat Shat;Shat' Rhat]
     exec LQR_problem_state_space.sci;
     

     //Passing arguments to the function to get as output the feedback matrices
     [Fphat,Fdhat,ns,nf,expect_cl_poles] = LQR_problem_state_space(A,B,[Qhat Shat;Shat' Rhat]);
     
     //Computation of the state-feedback gain matrices in the original basis
     
     if ns + nf ~= 0 then
        //Construction of feedback gain matrix for proportional state-feedback
         Fp = F + [Fphat zeros(no_inputs,dim_fast_system)]*inv(T);
         //Construction of feedback gain matrix for derivative state-feedback
         Fd = [Fdhat zeros(no_inputs,dim_fast_system)]*inv(T);
     else
         Fp = [];
         Fd = [];
     end
endfunction

//Function to block diagonilize a square matrix using RRE form
function [P,T] = block_Ematrix(E)
    //Size of the square matrix E
    n = size(E,1);
    //Pre-multiplying matrix P that will take E to its RRE form
    P =rref([E,eye(n,n)])(:,n+1:2*n); 
    //RRE form of E
    E1 = P*E;
    ////Pre-multiplying matrix T that will take (P*E)' to its RRE form
    T =rref([E1',eye(n,n)])(:,n+1:2*n);
    //Post-multiplying matrix T that will take PE to a block diagonal form with identity matrix as the block matrix
    T = T';
endfunction 
