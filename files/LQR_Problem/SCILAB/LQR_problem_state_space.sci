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
// Description: This Scilab program finds feedback gain matrices to solve the LQR problem for state-space LTI systems

// Input:
//   - LQR cost matrix represented as [Q S; S' R] 
//   - System matrices A, B, where the state-space LTI system is given by (dx/dt) = Ax + Bu

// Output:
//   - The state-feedback gain matrices Fp and Fd such that the control law u = Fp x + Fd (dx/dt) solves the LQR problem
//   - The closed loop system matrices of the system after u = Fp x + Fd (dx/dt) is employed. 
//   - The closed loop system is Ecl (dx/dt) = Acl x 

// Assumptions to be true:
//   - Cost matrix should be positive semidefinite where the input cost matrix R may or may not be singular
//   - The system must be (A,B) stabilizable
//*********************************************************************************************************************//
function [Fphat,Fdhat,ns,nf,expect_cl_poles] = LQR_problem_state_space(A,B,LQR_cost_matrix)
    //Check for (A,B) stabilizability
    exec pbhTest.sci;
    is_stabilizable = pbhTest(A,B);
    
    if is_stabilizable==0 then
        disp('System is stabilizable');
        disp('Finding the controller gain matrices');
    
        //Computation of number of states (no_states) and inputs (no_inputs)
        [no_states,no_inputs] = size(B);
    
        //Partition the cost matrix LQR_cost_matrix = [Qhat Shat;Shat' Rhat]
        Qhat = LQR_cost_matrix(1:no_states,1:no_states);
        Rhat = LQR_cost_matrix(no_states+1:no_states+no_inputs,no_states+1:no_states+no_inputs);
        Shat = LQR_cost_matrix(1:no_states,no_states+1:no_states+no_inputs);
    
        //Changing the basis of the input space for the problem
    
        //Diagonalize Rhat to Rtilde = [0 0;0 Re](for singular R) and Rtilde = Re (for non-singular R)
        rank_input_cost_mat = rank(Rhat);
         
        [U,Rtilde] = spec(Rhat);
        
        //Extract the non-singular diagonal block of Rtilde 
        Re = Rtilde(no_inputs - rank_input_cost_mat + 1:no_inputs,no_inputs - rank_input_cost_mat + 1:no_inputs);
        
        //Partition B conforming to the partition of Rhat based on rank
        B1 = clean((B*U)(:,1:no_inputs-rank_input_cost_mat)); 
        B2 = clean((B*U)(:,no_inputs-rank_input_cost_mat + 1:no_inputs));
        
        //Partition Shat conforming to the partition of Rhat based on rank
        S1 = clean((Shat*U)(:,1:no_inputs-rank_input_cost_mat));
        S2 = clean((Shat*U)(:,no_inputs-rank_input_cost_mat + 1:no_inputs));
        
        if rank_input_cost_mat~=0 then
            //Construction of the various matrices that form the reduced Hamiltonian matrix pairs
            //Er = [I 0;0 0] and Hr = [Ar -L Br;-Qr -Ar' 0;0 Br' 0];
            Ar = clean(A - B2*inv(Re)*S2');
            Qr = clean(Qhat - S2*inv(Re)*S2');
            L = clean(B2*inv(Re)*B2');
            Br = B1;    
        else
            Ar = clean(A);
            Qr = clean(Qhat);
            L =  zeros(no_states,no_states);
            Br = B1;
        end 
            
        //Construction reduced Hamiltonian matrix pairs (Er,Hr)
        d = no_inputs - rank_input_cost_mat;
        Hr = [Ar -L Br;-Qr -Ar' zeros(no_states,d);zeros(d,no_states) Br' zeros(d,d)];
        Er = blockdiag(eye(2*no_states,2*no_states),zeros(d,d));
         
        
        //Computation of the dimension of good slow space and fast space
        
        //Initialization of polynomial
        s = poly(0,'s');
        //Dimension of good slow space
        ns = (degree(det(s*Er-Hr)))/2;
        //Dimension of fast space
        nf = no_states - ns;
        
        //Executing function to find a basis of the good slow supspace
        exec basis_for_good_slow_space.sci;
            
        //Calling function to find a basis of the good slow subspace 
        //Output: - Columns of Vlambda is a basis
        //        - Gamma is a diagonal matrix containing the Lambda-set elements. Further Er*Vlambda*Gamma - Hr*Vlambda
        [Vlambda,expect_cl_poles] = basis_for_good_slow_space(Er,Hr,ns); 
        // V1lambda,V2lambda,V3lambda are partitioned conforming to a partion of states, costates, and input
        V1lambda = Vlambda(1:no_states,:);
        V2lambda = Vlambda(no_states+1:2*no_states,:);
        V3lambda = Vlambda(2*no_states+1:size(Er,1),:);

        //Executing function to find a basis of nf dimensional fast supspace
        exec basis_for_fast_space.sci;
        
        //Construction of PD state-feedback gain matrices 
        
        if nf ~= 0 then
            //This part is used if the underlying LQR problem is singular
            
            //Calling function to find a basis of nf dimensional fast subspace 
            //Output: Columns of W is a basis
            W = basis_for_fast_space(Ar,Br,Qhat,S1,S2,Re,nf)

            //Construction of matrices to find maximal rank minimizing solution of LQR LMI 
            X1 = [V1lambda W];
            X2 = [V1lambda Br Ar*W(:,1:nf-d)];
            
            //Construction of maximal rank minimizing solution of LQR LMI
            Kmax = [V2lambda zeros(no_states,no_states-ns)]*inv(X1);
            
            //Construction of non-unique part of feedback gain matrices
            G = find_gain_matrices(Ar,L,Kmax,Br,V3lambda,X1,nf,d)
            
            //Construction of feedback gain matrix for proportional state-feedback
            Fphat = U*[[V3lambda G]*inv(X1);-inv(Re)*(S2' + B2'*Kmax)];
            
            //Construction of feedback gain matrix for derivative state-feedback
            Fdhat = U*[[zeros(d,ns) eye(d,d) -G(:,1:nf-d)]*inv(X2);zeros(size(U,2)-d,ns+nf)];
        else 
            //This part is used if the underlying LQR problem is regular
            
            //Construction of matrices to find maximal rank minimizing solution of LQR LMI 
            Kmax = [V2lambda]*inv(V1lambda);
            
            //Construction of feedback gain matrix for proportional state-feedback
            Fphat = U*[-inv(Re)*(S2' + B2'*Kmax)];
            
            //For the regular LQR problem - Derivative control is not require (Fd = 0)
            Fdhat = zeros(no_inputs,no_states);
        end 
    else    
        disp('The system is not (A,B) stabilizable - Problem is not solvable');
        expect_cl_poles = [];
        Fphat = [];
        Fdhat = [];
        ns = 0;nf = 0;
    end    
endfunction


function [G] = find_gain_matrices(Ar,L,Kmax,Br,V3lambda,X1,nf,d)
    G = zeros(d,nf)
    // Create the matrix T
    T = Ar - L * Kmax + Br * [V3lambda, G];

    // Check if T is non-singular
    while det(T) == 0
        // If T is singular, update G to make T non-singular
        G = rand(d, nf);  // You can choose a different strategy to update G

        // Recreate the matrix T with the updated G
        T = Ar - L * Kmax + Br * [V3lambda, G];
     end 
endfunction
