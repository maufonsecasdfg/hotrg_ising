#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import scipy as sc
import scipy.linalg as scl
import matplotlib.pyplot as plt
import scipy.linalg


def tensorsvd(T_input,leftlegs,rightlegs,D='infinity'):
    '''Reshapes a tensor T_input into a matrix with first index corresponding
    to leftlegs and second index corresponding to rightlegs. Takes SVD
    and outputs U, s, V. U and V are reshaped to tensors leftlegs x D and 
    D x rightlegs respectively.
    '''
    T = np.transpose(T_input,leftlegs+rightlegs)
    xsize = 1
    leftsize_register = []
    for i in range(len(leftlegs)):
        xsize *= T.shape[i]
        leftsize_register.append(T.shape[i])
    ysize = 1
    rightsize_register = []
    for i in range(len(leftlegs),len(leftlegs)+len(rightlegs)):
        ysize *= T.shape[i]
        rightsize_register.append(T.shape[i])
    T = np.reshape(T,(xsize,ysize))
    
    U, s, V = np.linalg.svd(T,full_matrices = False)
    
    if D != 'infinity' and D < len(s):
        s = np.diag(s[:D])
        U = U[:,:D]
        V = V[:D,:]
    else:
        D = len(s)
        s = np.diag(s)
        
    U = np.reshape(U,leftsize_register+[D])
    V = np.reshape(V,[D]+rightsize_register)
        
        
    return U, s, V

def tensoreig(T_input,leftlegs,rightlegs,D='infinity'):

    T = np.transpose(T_input,leftlegs+rightlegs)
    xsize = 1
    leftsize_register = []
    for i in range(len(leftlegs)):
        xsize *= T.shape[i]
        leftsize_register.append(T.shape[i])
    ysize = 1
    rightsize_register = []
    for i in range(len(leftlegs),len(leftlegs)+len(rightlegs)):
        ysize *= T.shape[i]
        rightsize_register.append(T.shape[i])
    T = np.reshape(T,(xsize,ysize))
	
    TT = np.matmul(T,np.transpose(T))
    
    if D != 'infinity':
        if D > xsize:
            D = xsize
        lambd, U = scipy.linalg.eigh(TT,eigvals=(xsize-D,xsize-1))
    else:
        lambd, U = scipy.linalg.eigh(TT)
        
    lambd = lambd[::-1]
    lambd = np.diag(lambd)
        
    U = np.reshape(U,leftsize_register+[D])        
        
    return U, lambd

def ncon(tensorlist,labellist):
    '''Given a number of tensors connected in a diagram, specify the legs of 
    each tensor with a label and output the contraction of the diagram.
    The labeling is done as 1,2,3,... for the connected legs that are going
    to be contracted and in the order of contraction. The open legs are 
    labeled as -1, -2, -3, ..., representing the order in which they will be 
    represented in the output tensor.
    '''

    maxlabel = max(max(labellist)) 
    label = 1

    while label <= maxlabel:
        tensor1_found = False
        tensor2_found = False
        #Find the tensors that have the label for the leg that is going to be contracted next
        for ten_ind in range(len(labellist)):
            for i in range(len(labellist[ten_ind])):
                if labellist[ten_ind][i] == label:
                    #if a leg connect a tensor with itself, tensor1 and tensor2 are the same
                    if labellist[ten_ind].count(label) == 2:
                        labels = [label]
                        tensor1 = tensorlist[ten_ind]
                        tensor1_number = ten_ind
                        tensor1_index_con = [i]
                        tensor1_found = True
                        for j in range(i+1,len(labellist[ten_ind])):
                            if labellist[ten_ind][j] == label:
                                tensor2 = tensorlist[ten_ind]
                                tensor2_number = ten_ind
                                tensor2_index_con = [j]
                                tensor2_found = True
                                break
                        break
                    #if the labels connect different tensors, tensor1 and tensor2 are not the same
                    else:
                        if not tensor1_found:
                            tensor1 = tensorlist[ten_ind]
                            tensor1_number = ten_ind
                            tensor1_index_con = [i]
                            tensor1_found = True
                        else:
                            tensor2 = tensorlist[ten_ind]
                            tensor2_number = ten_ind
                            tensor2_index_con = [i]
                            tensor2_found = True
                        labels = [label]
                        #look if next label to contract is in the two chosen tensors, if it is, add it to the labels to contract
                        #keep looking until the condition is not met
                        look_for_more_contrac = True
                        while look_for_more_contrac and tensor2_found:
                            if label+1 in labellist[tensor1_number] and label+1 in labellist[tensor2_number]:
                                tensor1_index_con.append(labellist[tensor1_number].index(label+1))
                                tensor2_index_con.append(labellist[tensor2_number].index(label+1))
                                labels.append(label+1)
                                label += 1
                            else:
                                look_for_more_contrac = False
                        if tensor2_found:
                            break
            if tensor2_found:
                break
        
        #contract the two tensors in the indexes found if they are not the same tensor
        if tensor1_number != tensor2_number:
            tensor_contraction = np.tensordot(tensor1,tensor2,axes=(tensor1_index_con,tensor2_index_con))
            #update tensorlist so that the two just contracted tensors do not appear anymore
            tensorlist = tensorlist[:tensor1_number] + tensorlist[tensor1_number+1:tensor2_number] + tensorlist[tensor2_number+1:]
            
        #perfom a trace over the given indices if the labels belonged to the same tensor
        else:
            tensor_contraction = np.trace(tensor1,axis1=tensor1_index_con[0],axis2=tensor2_index_con[0])
            #update tensorlist so that the tensor does not appear anymore (traced version to be added back soon)
            tensorlist = tensorlist[:tensor1_number] + tensorlist[tensor1_number+1:]
            
        #construct the leg labels for the new tensor resulting from the contraction
        A = []
        for lab in labellist[tensor1_number]:
            if lab not in labels:
                A.append(lab)        
        labellist[tensor1_number] = A[:]
        if tensor1_number != tensor2_number:
            B = []
            for lab in labellist[tensor2_number]:
                if lab not in labels:
                    B.append(lab) 
            labellist[tensor2_number] = B[:]
        
            new_tensor_labels = labellist[tensor1_number] + labellist[tensor2_number]            
            labellist = labellist[:tensor1_number] + labellist[tensor1_number+1:tensor2_number] + labellist[tensor2_number+1:]

            
        else:
            new_tensor_labels = labellist[tensor1_number]
            labellist = labellist[:tensor1_number] + labellist[tensor1_number+1:]

        #add the resulting tensor and its labels to the begining of tensorlist and labellist
        tensorlist.insert(0, tensor_contraction)
        labellist.insert(0, new_tensor_labels)
        
        label += 1
    
    contracted_tensor = tensorlist[0]
    indices = labellist[0]
    
    #order the tensor indices so that they appear according to the -1,-2-3,... order
    indices = [-i for i in indices]
    
    open_leg_order = []
    for i in range(1,len(indices)+1):
        for j in range(len(indices)):
            if i == indices[j]:
                open_leg_order.append(j)
    
    contracted_tensor = np.transpose(contracted_tensor,open_leg_order)
        
    return contracted_tensor
        
        

