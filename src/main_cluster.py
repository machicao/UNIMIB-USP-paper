#!/usr/bin/python3
 
import networkx as nx
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join

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
            #atualizacao  dos arrays
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
     
def extract_measures(G):
    """Extract a vector of measurements
    """
    num_nodes =  nx.number_of_nodes(G)
    if num_nodes == 0:
        raise Exception('graph is empty')   
    else: 
        F = {} #features 

        """ Global """
        H = hierarchical_degree(G, nivel_max=3)
        F['avg_degree'] =  np.average(H[0,:])
        F['avg_hier2']  =  np.average(H[1,:])
        F['avg_hier3']  =  np.average(H[2,:])
        #F['avg_clust']  =  nx.average_clustering(G) #0, para este caso todas as redes obtiver 0
        #tic()  #avgpathlenght is expensive
        F['avg_path']   =  nx.algorithms.average_shortest_path_length(G)
        F['pearson']    =  nx.algorithms.degree_pearson_correlation_coefficient(G) #-0.18656860948991077   
        #F['avg_leff']   =  nx.local_efficiency(G) 
        #F['avg_geff']   =  nx.global_efficiency(G)
        #F['avg_trans']  =  nx.algorithms.transitivity(G) #0

        """ Local measures """
        #ecc  =   nx.algorithms.eccentricity(G)          # [0.01408451 0.   0.32394366 0.   0.66197183]
        #tic() #betweeness is expensive
        #bet  =   nx.algorithms.betweenness_centrality(G)# [0.97183099 0.   0.   0.01408451 0.01408451]
        #katz =   nx.algorithms.katz_centrality(G)       # [0.67605634 0.29577465 0.    0.  0.02816901]
        #clust  = nx.algorithms.clustering(G)           #[0. 0. 1. 0. 0.] O CLUSTERING NESTAS REDES DA O
        #sqclust= nx.algorithms.square_clustering(G)    #[0.64788732 0.04225352 0.02816901 0.01408451 0.26760563]
        #close  = nx.algorithms.closeness_centrality(G) #[0.64788732 0.32394366 0.    0.    0.02816901]  
        #degcen = nx.algorithms.degree_centrality(G)    #[0.97183099 0.   0.     0.    0.02816901]


        #F['avg_ecc']    = np.average(list(ecc.values()))
        #F['avg_bet']    = np.average(list(bet.values()))    
        #F['avg_katz']   = np.average(list(katz.values()))
        #F['avg_sqclust']= np.average(list(sqclust.values()))
        #F['avg_close']  = np.average(list(close.values()))
        #F['avg_degcen'] = np.average(list(degcen.values()))

        return F 

# def histo_descriptor(array, bins):
#     hist, bin_edges = np.histogram(array, bins = bins)
#     return hist/sum(hist) #density_histo


#reactions = pd.read_csv('UNIMIB_data/ReactionMetabolites_list/R22_reactions.txt',sep="\t")
metabolites= pd.read_csv('UNIMIB_data/ReactionMetabolites_list/R22_metabolites.txt',sep="\t")


def read_network(filename): 
    """ return complete gcc network"""
    df = pd.read_csv(filename,sep="\t")
    nodeA= df.iloc[:,0]
    nodeB= df.iloc[:,1]
    peso = df.iloc[:,2]

    all_nodes =pd.concat([nodeA, nodeB], ignore_index =True) #ignore_index, i.e re-indexa

    G=nx.Graph() #DiGraph

    #add nodes
    #for m in all_nodes[all_nodes.isin(metabolites['ID'])]:
    #    #G.add_node(m, color='r') #metabolite    
    #    G.add_node(m) #metabolite    
    #for r in all_nodes[~all_nodes.isin(metabolites['ID'])]:
    #    #G.add_node(r, color='b') #reaction
    #    G.add_node(r) #reaction
    #add edges
    for i in range(nodeA.size):
        e1, e2, p=nodeA.get(i), nodeB.get(i), peso.get(i)
        G.add_edge(e1, e2, peso=float(p))   
    return G    
    # largest connected component
    #G0 = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]#0 = the largest network
    #return G0

def threshold_network_range(G, thre_min, thre_max):
    """ filtering by threshold """
    filtered=[(u,v) for (u,v,w) in G.edges(data=True) if (w['peso'] >=thre_min and w['peso'] < thre_max)]

    if(len(list(filtered)) == 0):
        raise Exception('filtered graph is empty')  
    else: 
        # largest connected component
        Gthre = sorted(nx.connected_component_subgraphs(nx.Graph(filtered)), key=len, reverse=True)[0]
        #print('thre= [%.4f, %.4f], nodes= %d'% (thre_min, thre_max, nx.number_of_nodes(Gthre)))
        return Gthre        
    #write
    # nx.write_graphml(G0_thre, 'TCGA_A7_A0CE_thre=0.7_0.8.graphml')  
    
def threshold_network_less_than(G, thre):
    """ filtering by threshold """
    filtered=[(u,v) for (u,v,w) in G.edges(data=True) if (w['peso'] <= thre)]
    if(len(list(filtered)) == 0):
        raise Exception('filtered graph is empty')  
    else: 
        # largest connected component
        Gthre = sorted(nx.connected_component_subgraphs(nx.Graph(filtered)), key=len, reverse=True)[0]
        #print('thre= [%.4f, %.4f], nodes= %d'% (thre_min, thre_max, nx.number_of_nodes(Gthre)))
        return Gthre         
    
def threshold_network_bigger_than(G, thre):
    """ filtering by threshold """
    filtered=[(u,v) for (u,v,w) in G.edges(data=True) if (w['peso'] >= thre)]

    if(len(list(filtered)) == 0):
        raise Exception('filtered graph is empty')  
    else: 
        # largest connected component
        Gthre = sorted(nx.connected_component_subgraphs(nx.Graph(filtered)), key=len, reverse=True)[0]
        #print('thre= [%.4f, %.4f], nodes= %d'% (thre_min, thre_max, nx.number_of_nodes(Gthre)))
        return Gthre        
    #write
    # nx.write_graphml(G0_thre, 'TCGA_A7_A0CE_thre=0.7_0.8.graphml')      


def main(labels):
	# Load network
	#left-threshold
	# thre= [0.0, 10**-4, 10**-3, 10**-2, 0.1, 0.2,  0.3 ,0.4, 0.5,0.6, 1.1]   
	# thre= [0.0, 0.1, 0.2,  0.3 ,0.4, 0.5,0.6, 1.1]   
	# thre= [0.0, 10**-4, 10**-3]   

	#bigger_than threshold
	#thre= [0.0, 10**-4, 10**-3, 10**-2, 0.1, 0.2,  0.3 ,0.4, 0.5,0.6, 0.7, 0.8, 0.9]
	thre= [10**-4, 10**-3,10**-2, 0.1, 0.2,  0.3 ,0.4, 0.5,0.6, 0.7]

	pacientes= [f for f in listdir('approaches/weigth/cancer') if isfile(join('approaches/weigth/cancer', f))] 
	#cancer has the same pacientes as normal

	path= 'approaches/weigth/'
	#labels= ['cancer', 'normal']

	num_pacientes=len(pacientes)
	num_medidas=5 #'avg_degree', 'avg_hier2','avg_hier3','avg_path','pearson'
	#num_thre=(len(thre)-1) #Range threshold
	num_thre =len(thre)#>= threshold   
	num_descritores= num_thre * num_medidas

	features =  np.zeros((num_pacientes*2, num_descritores))
	k = list()
	str_clase=list()
	int_clase=list()
	atr_names= ['avg_degree', 'avg_hier2','avg_hier3','avg_path','pearson']*num_thre
	for c in range(len(labels)):
	    print(labels[c])
	    for p in range(num_pacientes):
	        G= read_network(path+labels[c]+'/'+ pacientes[p]) #reading 
	        print(pacientes[p])
	        k.append(pacientes[p])
	        str_clase.append(labels[c])
	        int_clase.append(c)
	        #threshold      
	        #v=len(thre)-1 # #Range threshold
	        v= len(thre)  #>= threshold   
	        for i in range(v):
	            #Gthre=threshold_network_range(G, thre[i], thre[i+1]) # #Range threshold
	            Gthre=threshold_network_bigger_than(G, thre[i]) #>= threshold     
	            num_nodes =  nx.number_of_nodes(Gthre)
	            print(num_nodes)
	            f=extract_measures(Gthre)  #dict 
	            feats= list(f.values())
	            #atr_names= list(f.keys())
	                        
	            features[(c*num_pacientes)+p][i*num_medidas:(i*num_medidas)+num_medidas] = feats[:]             
	            #print('%s, thre= [%.4f, %.4f], %s'% (pacientes[p].replace('.txt',''), thre[i], thre[i+1],feats[:] ))
	#salvar pd            
	atr_names= ['avg_degree', 'avg_hier2','avg_hier3','avg_path','pearson']*num_thre
	# print(atr_names)
	data = pd.DataFrame(features, columns=atr_names)
	col1 = pd.DataFrame(k,columns= ['paciente'])
	col2 = pd.DataFrame(str_clase,columns= ['name_class'])
	col3 = pd.DataFrame(int_clase,columns= ['class'])

	df = pd.concat([data, col1['paciente'], col2['name_class'],col3['class']], axis = 1)
	df = pd.concat([data,   col3['class'] ], axis = 1)

	# df

	#save dataframe to csv 

	df = pd.concat([data, col1['paciente'], col2['name_class'],col3['class']], axis = 1)
	df.to_csv('biggerthan'+ labels[0] +'.csv', sep='\t', index=False)


if __name__ == '__main__':    
#   labels= ['cancer'] #PARTE1
    labels= ['normal'] #PARTE2
    main(labels)
