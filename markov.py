# Markov's Clustering Algorithm

import numpy as np
from numpy import linalg as la

# Loads the corresponding file into the association matrix
def load_dataset(file_name):
    
    # Loading the ATT file.
    if(file_name == "attweb_net.txt"):
        data = np.loadtxt(file_name, delimiter=' ', dtype = int)        
        myList = []
        for i in data:
            if i[0] not in myList:
                myList.append(i[0])
            if i[1] not in myList:
                myList.append(i[1])
        n = len(myList)
        s = (n,n)
        association_matrix = np.zeros(s)
        
        #Populate the Association Matrix
        for i in data:
            p = myList.index(i[0])
            q = myList.index(i[1])
            association_matrix[p][q] = 1
            association_matrix[q][p] = 1
        
        # Adding self loops
        association_matrix[range(n), range(n)] = 1
        
    # Loading the Physics file.
    if(file_name == "physics_collaboration_net.txt"):
        data = np.loadtxt(file_name, delimiter=' ', dtype = 'str')
        myList = []
        for i in data:
            if i[0] not in myList:
                myList.append(i[0])
            if i[1] not in myList:
                myList.append(i[1])
        n = len(myList)
        s = (n,n)
        association_matrix = np.zeros(s)
        
        #Populate the Association Matrix
        for i in data:
            p = myList.index(i[0])
            q = myList.index(i[1])
            association_matrix[p][q] = 1
            association_matrix[q][p] = 1
            
        # Adding self loops
        association_matrix[range(n), range(n)] = 1
    
    # Loading the Yeast file.
    if(file_name == "yeast_undirected_metabolic.txt"):
        data = np.loadtxt(file_name, delimiter="\t", dtype = int)
        myList = []
        for i in data:
            if i[0] not in myList:
                myList.append(i[0])
            if i[1] not in myList:
                myList.append(i[1])
        n = len(myList)
        s = (n,n)
        association_matrix = np.zeros(s)
        
        #Populate the Association Matrix
        for i in data:
            p = myList.index(i[0])
            q = myList.index(i[1])
            association_matrix[p][q] = 1
            association_matrix[q][p] = 1
            
        # Adding self loops
        association_matrix[range(n), range(n)] = 1
    return association_matrix, n
    
# Performs normalization on the given association matrix
def normalize_association_matrix(association_matrix):
    x = association_matrix.sum(axis=0)
    association_matrix = association_matrix / x
    return association_matrix

# Expands the association matrix as per the specified power parameter 'e'
def expand_association_matrix(association_matrix, expansion_power):
    association_matrix = la.matrix_power(association_matrix, expansion_power)
    return association_matrix

# Inflates the association matrix as per the specified inflation parameter 'p'
def inflate_association_matrix(association_matrix, inflation_power):
    association_matrix = np.power(association_matrix, inflation_power)
    return association_matrix

# Sets the entries in the association matrix to zero whose value is smaller than the threshold value.
def round_off_values(association_matrix, threshold):
    for rows in range(association_matrix.shape[0]):
        for cols in range(association_matrix.shape[1]):
            x = association_matrix[rows][cols]
            if(x <= threshold):
                association_matrix[rows][cols] = 0
    return association_matrix

# Interprets cluster nodes from the resultant matrix
def form_clusters(resultant_matrix):
    cluster_nodes = set([])
    for rows in range(resultant_matrix.shape[0]):
        for cols in range(resultant_matrix.shape[1]):
            if(rows == cols):
                cluster_list = list()
                for position, value in enumerate(resultant_matrix[rows,:]):
                    if(value > 0):
                        node = position
                        cluster_list.append(node)
                my_tuple = tuple(cluster_list)
                cluster_nodes.add(my_tuple)
    return cluster_nodes
                
# Implements Markov's Clustering Algorithm on the association matrix
def markov_clustering(association_matrix, expansion_power, inflation_power):
    association_matrix = normalize_association_matrix(association_matrix)
    # specifies round off pruning threshold so as to set values less than threshold value directly to zero
    threshold = 0.005
    count = 0
    # specifies the maximum number of iterations - in case the MCL does not converge
    maximum_counter = 100
    previous_matrix = np.array(association_matrix)
    flag = True
    while(flag):
        association_matrix = expand_association_matrix(association_matrix, expansion_power)
        association_matrix = inflate_association_matrix(association_matrix, inflation_power)
        association_matrix = normalize_association_matrix(association_matrix)
        association_matrix = round_off_values(association_matrix, threshold)
        count = count + 1
        if((count > maximum_counter) or (np.array_equal(association_matrix, previous_matrix))):    
            print("Number of Iterations : " + str(count))
            flag = False
        previous_matrix = association_matrix
    return association_matrix

# Writes the cluster nodes to the corresponding .clu file which acts as an input for visualization using PAJEK.
def write_cluster_file(file_name, cluster_nodes, vertices):
    splitter = file_name.split("_")
    output_file_name = splitter[0] + ".clu"
    f = open(output_file_name,"w") 
    f.write("*Vertices " + str(vertices) + "\n")
    print("Number of Vertices : " + str(vertices))
    cluster_hash_map = {}
    print("Number of Clusters : " + str(len(cluster_nodes)))
    for position, value in enumerate(cluster_nodes):
        for i, j in enumerate(value):
            pos = position + 1
            cluster_hash_map[j] = pos
    
    for key, value in cluster_hash_map.items():
        f.write(str(value) + "\n")
    f.close()

    
# Specify the File name(uncomment) on which you want to run Markov Clustering Algorithm and comment the other two files
file_name = "attweb_net.txt"
#file_name = "physics_collaboration_net.txt"
#file_name = "yeast_undirected_metabolic.txt"

# Load the corresponding dataset in order to form association matrix
association_matrix, vertices = load_dataset(file_name)

# Specify the expansion / power parameter 'e'
expansion_parameter = 4

# Specify the inflation parameter 'r'
inflation_parameter = 1.8

# Invoke the Markov Clustering function on the specified parameters
association_matrix = markov_clustering(association_matrix, expansion_parameter, inflation_parameter)

# Interpret clusters from the resultant matrix
cluster_nodes = form_clusters(association_matrix)

# Write the cluster nodes into the corresponding .clu file which acts as an input for PAJEK 
write_cluster_file(file_name, cluster_nodes, vertices)
