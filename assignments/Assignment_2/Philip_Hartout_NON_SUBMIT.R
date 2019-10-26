#########################################################
#######        COMPUTATIONAL BIOLOGY         ############
#######             HOMEWORK 2               ############
#########################################################
#                                                       #
# Simulate the evolution of sequences on a fixed tree   #
# under the TN93 sequence evolution model               #
#                                                       #
#########################################################
#########################################################

#########################################################
######    Code below this should not be changed   #######
#########################################################

library(Matrix)
library(ape)

nucleotides = c("T", "C", "A", "G")

get_nucleotide_from_number = function(i) {
    # Transform a nucleotide number into the appropriate letter
    return(nucleotides[i])
}

transform_to_nucleotides = function(sequence) {
    # Transform a sequence of nucleotide numbers into the appropriate character sequence
    nucl_sequence = paste0(lapply(sequence, get_nucleotide_from_number), collapse = "")
    return(nucl_sequence)
}

#########################################################
######    Code above this should not be changed   #######
#########################################################

# In all functions the following parameters are the same:
# pi: the stationary frequencies of nucleotides
# alpha1: rate coefficient for the C <-> T transition
# alpha2: rate coefficient for the A <-> G transition
# beta: rate coefficient for transversions
# N: number of sites in the simulated alignment

create_TN93_Q_matrix = function(pi, alpha1, alpha2, beta) {
    # Create the TN93 transition rate matrix Q as specified in the assignment.
    
    Q = matrix(c(-(alpha1*pi[2] + beta*pi[3] + beta*pi[4]), alpha1*pi[2], beta*pi[3], beta*pi[4],
      alpha1*pi[1], -(alpha1*pi[1] + beta*pi[3] + beta*pi[4]), beta*pi[3], beta*pi[4],
      beta*pi[1], beta*pi[2], -(beta*pi[1] + beta*pi[2] + alpha2*pi[4]), alpha2*pi[4],
      beta*pi[1], beta*pi[2], alpha2*pi[3], -(beta*pi[1] + beta*pi[2] + alpha2*pi[3])), 
      nrow=4,              # number of rows 
      ncol=4,              # number of columns 
      byrow = TRUE)       # fill matrix by rows 
    
    # Return the transition rate matrix
    # Q: 4 by 4 matrix of rates
    return(Q)
}

get_starting_nucleotide = function(pi) {
    # Sample a starting nucleotide from the stationary distribution
    
    nucleotide = sample(c(1,2,3,4), size = 1, prob = pi)

    # Return the sampled nucleotide
    # nucleotide: integer nucleotide value
    return(nucleotide)
}

get_starting_sequence = function(pi, N) {
    # Sample a starting sequence of length N
    starting_sequence = vector(mode="integer", length=N)
    for (i in 1:N) {
      starting_sequence[i] = get_starting_nucleotide(pi)
    }
    # Return the sampled sequence
    # starting_sequence: vector of integer nucleotide values
    return(starting_sequence)
}

get_evolved_sequence = function(sequence, branch_length, Q) {
    # Evolve a given nucleotide sequence along a branch of specified length.
    # sequence: nucleotide sequence at the beginning of the branch
    # branch_length: the length of the branch along which evolution happens
    # Q: the transition rate matrix
    
    probability_matrix = expm(branch_length*Q)
    print(probability_matrix)
    evolved_sequence = vector(mode="integer", length=length(sequence))
    
    for (i in 1:length(sequence)) {
      evolved_sequence[i] = sample(c(1,2,3,4), size = 1, prob = probability_matrix[sequence[i],])
    }
    
    # Return the nucleotide sequence after all positions have evolved along the given branch.
    # evolved_sequence: the vector of new integer nucleotide values at the end of the branch
    return(evolved_sequence)
}


simulate_evolution = function(newick_tree, pi, alpha1, alpha2, beta, N) {
    # Simulate evolution along the given tree.
    # newick_tree: the tree in newick text format
    
    # Transfrom the tree from text format to an object of the phylo class which represents the tree in R
    tree = read.tree(text = newick_tree)
    # Reorder the tree for easier traversing
    tree = reorder(tree, order = "cladewise")
    
    # Set up the Q matrix
    Q = create_TN93_Q_matrix(pi, alpha1, alpha2, beta)
    
    # Set up the starting sequence @ the root of the tree
    starting_sequence = get_starting_sequence(pi, N)
    
    # Prepare a list to store evolved sequences at each node
    sequence_per_node = list()
    sequence_per_node[[tree$edge[1,1]]] = starting_sequence
    
    # Walk the tree while evolving sequences along appropriate branches
    for (i in 1:length(tree$edge.length)) {
        node_parent = tree$edge[i, 1]
        node_child = tree$edge[i, 2]
        branch_length = tree$edge.length[i]
        parent_sequence = sequence_per_node[[node_parent]]
        
        child_sequence = get_evolved_sequence(parent_sequence, branch_length, Q)

        sequence_per_node[[node_child]] = child_sequence
    }

    # Transform the alignment from nucleotide indices to nucleotide characters
    # and filter out the sequences at the tips
    alignment = list()
    for (i in 1:length(tree$tip.label)) {
        alignment[[tree$tip.label[i]]] = transform_to_nucleotides(sequence_per_node[[i]])
    }
    
    # Return the simulated alignment.
    # The alignment should be in the form of a list where the tip label corresponds to the
    # appropriate simulated sequence, e.g. alignment$human = ACTG
    return(alignment)
}

test_simulation = function() {
  library(ape)
  newick_tree = "(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75);"
  N = 40
  beta = 0.035
  alpha1 = 0.044229
  alpha2 = 0.021781
  pi = c(0.22, 0.26, 0.33, 0.19)
  
  result = simulate_evolution(newick_tree, pi, alpha1, alpha2, beta, N)  
  return(result)
}
beta = 0.035
alpha1 = 0.044229
alpha2 = 0.021781

pi = c(0.22, 0.26, 0.33, 0.19)
sequence = get_starting_sequence(pi, 40)
Q = create_TN93_Q_matrix(pi, alpha1, alpha2, beta)

for (i in 1:9) {
  t = i * 100
  print(t)
  get_evolved_sequence(sequence, t, Q)
}

# library (stringr)
# emprical_pi = c(0,0,0,0)
# continue = F
  # x = 100 
# within_range = c(F, F, F, F)
# while (continue == F) {
#   result = test_simulation(x)
#   T_freq = str_count(result$human, "T")/nchar(result$human)
#   C_freq = str_count(result$human, "C")/nchar(result$human)
#   A_freq = str_count(result$human, "A")/nchar(result$human)
#   G_freq = str_count(result$human, "G")/nchar(result$human)
#   empirical_pi = c(T_freq, C_freq, A_freq, G_freq)
#   within_range_2 = c(F, F, F, F)
#   for (element in 1:length(emprical_pi)) {
#     upper_bound = pi[element]+10e-4
#     lower_bound = pi[element]-10e-4 
#     if (empirical_pi[element] > lower_bound & empirical_pi[element] < upper_bound) {
#       within_range[element] = T
#       within_range_2[element] = T
#     } 
#   }
#   x = x + 100
#   print(within_range_2)
#   print(x)
#   print(within_range)
#   if (all(within_range) == T){
#     continue = T
#   }
# }
# print(x)
