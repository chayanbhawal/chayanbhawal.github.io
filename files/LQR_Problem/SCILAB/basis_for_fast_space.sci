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
//  Description: This Scilab program finds basis vectors for a fast subspace of an output-nulling Hamiltonian system. 
//
//  Input:
//  - Partitioned cost-matrix [Qhat S1 S2;S1' 0 0;S2' 0 Re] where Re is positive definite
//  - Dimension of fast subspace: nf
//
//  Output:
//  - A matrix Vlambda whose columns are a basis of  good slow subspace of a output-nulling Hamiltonian system 
//**********************************************************************************************************************//
function W = basis_for_fast_space(Ar,Br,Qhat,S1,S2,Re,nf)
    //d = no_inputs - rank_input_cost_mat
    d = size(S1,2);
    
    //Rank of input cost matrix
    rank_input_cost_mat = rank(Re);
    
    //Cholesky factorization of [Qhat S1 S2;S1' 0 0;S2' 0 Re] to find C,D2 such that C'C = Qhat, D2'D2 = Re, C'D2 = S2
    
    //SVD of [Qhat S1 S2;S1' 0 0;S2' 0 Re] 
    [M,Z,N] = svd([Qhat S1 S2;S1' zeros(d,d) zeros(d,rank_input_cost_mat);S2' zeros(rank_input_cost_mat,d) Re],'econ');
    
    // Find T matrix such that [Q S1 S2;S1' 0 0;S2' 0 Re] = T'*T 
    T = sqrt(clean(Z))*N';                                              
    
    //Rank of [Qhat S1 S2;S1' 0 0;S2' 0 Re] 
    p = rank(clean(Z));
    
    //Partitioning T to construct C such that C'C = Qhat
    C = T(1:p,1:no_states);
    
    //Partitioning T to construct D2 such that D2'*D2 = Re and C'D2 = S2
    D2 = T(1:p,no_states+d+1:no_states+no_inputs);

    if rank_input_cost_mat~=0 then
        Cr = C - D2*inv(Re)*S2';
    else
        Cr = C;
    end

    //Executing function to form the Markov matrix required to find basis for fast subspace
    exec markov.sci;
    //Calling function to find the Markov matrix
    M = markov(Ar,Br,Cr,nf);
    N = kernel(M);
    Mat = [];
    if nf > d then
        for i = 1:nf-d
            Mat = [Mat Ar^(i-1)*Br];
        end
    else
        Mat = Br;
    end

    Wtilde = Mat*N;
    if nf > d then
        exec extend_basis.sci;
        W = extend_basis(Wtilde,[Br Ar*Wtilde]);
    else
        W = Wtilde;
    end
endfunction
