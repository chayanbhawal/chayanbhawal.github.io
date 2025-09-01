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

//*********************************************************************************************************************//
// Description: This Scilab program computes PD state-feedback closed-loop gain matrices to solve the LQR problem for both state-space and DAE systems. 
//
// Input:
//   - LQR cost matrix represented as [Q S; S' R]
//   - System matrices E, J, L, where the DAE system is given by E(dx/dt) = Jx + Lu
//     For the state-space system, E is the identity matrix.
//
// Output:
//   - Feedback gain matrices Fp, Fd where u = Fp * x + Fd * (dx/dt) that solves the LQR problem.
//   - Dimension of the good slow subspace (dim_good_slow) and fast subspace (dim_fast) of the (underlying) LQR problem.
//   - Expected closed loop poles (expect_cl_poles)
// 
// Assumptions:
//   - For DAE problem: The system must be impulse controllable.
//   - For standard state-space problem: The system must be stabilizable.
//**********************************************************************************************************************//

//Clear screen
clc
//Clear all variables
clear

//Executing function to load system and cost matrix
exec input_data_matrix.sci;                                             
//Load the system matrices 
//For DAE based LQR problem system matrix is E dx/dt = Jx + Lu
[E,J,L,LQR_cost_matrix] = input_data_matrix();
s = poly(0,'s');

//Check for regularity of the system
is_regular = 1;
if det(s*E-J)==0 then
    is_regular = 0;
    disp('The system is not regular - out of scope of this program.');
else
    disp('The system is regular - proceeding to check +ve semi-definiteness of cost matrix');
    //Check for positive semi-definiteness of cost matrix
    exec is_positive_semidefinite_matrix.sci;
    is_positive_semidefinite =  is_positive_semidefinite_matrix(LQR_cost_matrix);
    if is_positive_semidefinite == 1
        disp('The cost matrix is +ve semi-definite - proceeding to check type of the system');
        
        //Computation of number of states (no_states) and inputs (no_inputs)
        [no_states,no_inputs] = size(L);
    
        //Checking if the problem is DAE based
        if rank(E) < size(E,1) then
             flag = 0;
             disp('The system is a DAE - proceeding to check for the validity of the assumptions');
             
             //Executing function to solve LQR problem for DAE
             exec LQR_problem_DAE.sci; 
            
             //Passing arguments to the function to get as output the feedback matrices
             [Fp,Fd,dim_good_slow,dim_fast,expect_cl_poles] = LQR_problem_DAE(E,J,L,LQR_cost_matrix,no_states,no_inputs);
    
        else if E == eye(no_states,no_states)
             flag = 1
             disp('The system is in standard state-space form  - proceeding to check for the validity of the assumptions');
        
             //Executing function to solve LQR problem for state-space
             exec LQR_problem_state_space.sci;                                        
        
             //Passing arguments to the function to get as output the feedback matrices
             [Fp,Fd,dim_good_slow,dim_fast,expect_cl_poles] = LQR_problem_state_space(J,L,LQR_cost_matrix);
        
        else
            flag = 1;
            disp('The system is not a DAE - reframing the problem as a standard state-space problem and solving') 
        
            //Executing function to solve LQR problem for state-space
            exec LQR_problem_state_space.sci;                                        
        
            //Passing arguments to the function to get as output the feedback matrices
            [Fp,Fd,dim_good_slow,dim_fast,expect_cl_poles] = LQR_problem_state_space(inv(E)*J,inv(E)*L,LQR_cost_matrix);
        end
    end
        if dim_good_slow + dim_fast ~= 0
        //Construction of the closed loop system matrices Ecl (dx/dt) = Acl x 
         Ecl = E - L*Fd;
         Acl = J + L*Fp;
         
         cl_system = clean(roots(det(s*Ecl - Acl)));
        
        //Displaying all the results on the command console
        disp('Good slow subspace dimension',dim_good_slow);
        if expect_cl_poles == [] then
            disp('The closed loop system will not admit any slow subspace');
        else 
            disp('The expected closed loop poles are: ', expect_cl_poles');
        end
    
        disp('Fast subspace dimension',dim_fast);
        if flag == 0 & dim_fast == 0
            disp('The underlying state-space LQR problem is regular - P controller required to solve the LQR problem')
            disp('The control law to solve the problem is  u = Fp * x with');
            disp('Fp = ',clean(Fp)); 
        elseif flag == 1 & dim_fast == 0
            disp('The state-space LQR problem is regular - P controller required to solve the LQR problem')
            disp('The control law to solve the problem is  u = Fp * x with');
            disp('Fp = ',clean(Fp)); 
        elseif flag == 0 & dim_fast == 1
            disp('The underlying state-space LQR problem is singular - PD controller required to solve the LQR problem')
            disp('The control law to solve the problem is  u = Fp * x + Fd * (dx/dt) with');
            disp('Fp = ',clean(Fp),'Fd = ',clean(Fd)); 
        else flag == 1 & dim_fast == 1
            disp('The state-space LQR problem is singular - PD controller required to solve the LQR problem')
            disp('The control law to solve the problem is  u = Fp * x + Fd * (dx/dt) with');
            disp('Fp = ',clean(Fp),'Fd = ',clean(Fd)); 
        end
        
        disp("The closed loop system is Ecl(dx/dt) = Acl x with");
        disp("Ecl =",clean(Ecl),"Acl=",clean(Acl));
        
        if cl_system ~= [] then
            disp("The roots of the closed loop system on application of feedback",cl_system);
        else
            disp('The closed loop system is a DAE with only fast subsystem');
        end
        end
    else
            disp('The cost matrix is not +ve semi-definite - solution out of scope of this program');
    end
end
