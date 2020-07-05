import numpy as np
import networkx as nx

def hierarchical_degree(G, nivel_max):

    """Returns a feature vector (matrix) containing the 
    degree of each level
    Example of use: 
    res = hierarchical_degree(G, nivel_max=3)
    degree= res[0,:]
    hier2=res[1,:]
    hier3=res[2,:]
    """
    listNodes = G.nodes()
    num_nodes=len(list(listNodes))
    features =  np.zeros((nivel_max,num_nodes))
    #print(num_nodes) 
    for item ,n in zip(listNodes, range(0,num_nodes)):
        #print(n,'/',len(listNodes))
        #inicializacao
        vizinhos = nx.neighbors(G,item)
        lista_vizinhos=list(vizinhos) 
        root = [item] * len(lista_vizinhos)
        ray = [1] * len(lista_vizinhos) 
        #print(len(lista_vizinhos))
        while(len(lista_vizinhos)  > 0):  
            #get the element
            theVizi, theRoot, theRay = lista_vizinhos[0], root[0], ray[0] 
            features[theRay-1][n] = float(features[theRay-1][n]) + 1
            #print(features[theRay][n])         
            #atualização  dos arrays
            if(theRay < nivel_max):
                vizinhosVizinhos =  G.neighbors(theVizi)
                lista_vizinhosVizinhos= list(vizinhosVizinhos)
                vizinhosVizinhos = list(set(lista_vizinhosVizinhos) - set(lista_vizinhos))
                lista_vizinhos = lista_vizinhos + vizinhosVizinhos
                rootRoot = [theVizi] * len(vizinhosVizinhos)
                root = root + rootRoot
                rayRay = [theRay + 1] * len(vizinhosVizinhos) 
                ray = ray + rayRay
            del lista_vizinhos[0]
            del root[0]
            del ray[0]           
    return features
    #print(features) 