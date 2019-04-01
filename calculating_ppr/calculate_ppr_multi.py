from multiprocessing import Process, Queue
import numpy as np
import pandas as pd        
import pickle 
import numpy as np
import scipy.sparse
import igraph as ig


import argparse

parser = argparse.ArgumentParser(description='Make translate dict')
parser.add_argument('-l', '--language', type=str,
                    help='Laugnage code of language')
parser.add_argument('-a', '--alpha', type=float,
                    help='Alpha of PPR calculation')
parser.add_argument('-k', '--k_factor', type=float,
                    help='k factor that determine random-walker behavior')




def calculate_ppr(core_number, node_list, W_matrix, alpha, output):
    total_length = W_matrix.shape[0]
    result_dict = {}
    for count, node in enumerate(node_list):
        e_s = scipy.sparse.csr_matrix(([1], ([0],[node])), shape=(1, total_length))
        p = scipy.sparse.csr_matrix(([1], ([0],[node])), shape=(1, total_length))
        for i in range(100):
            new_p =  (alpha * e_s) + ((1 - alpha) * (p * W_matrix))
            diff = abs(new_p - p)
            if diff.max() < 0.0001:
                p = new_p
                break   
            p = new_p
        p =  normalize_and_convert_to_dict(p)
        result_dict[node] = p
        if (count + 1)%1000 == 0:
            print('core ' + str(core_number) + ': ' + str(count/len(node_list)*100) + ' % completed') 
    output.put(result_dict)        
    print('core ' + str(core_number) + ' process_completed')   

    
def normalize(d, target=1.0):
    raw = sum(d.values())
    factor = target/raw
    return {key:value*factor for key,value in d.items()}

def normalize_and_convert_to_dict(csr_matrix):
    dummy_dict = csr_matrix.todok()
    dummy_dict = normalize(dummy_dict)
    return {k[1]:v for k,v in dummy_dict.items() if v>0}
    
    

if __name__ == '__main__':
    args = parser.parse_args()
    lang_code = args.language
    alpha = args.alpha
    k_factor = args.k_factor
    number_of_process = 30
    
    G = ig.Graph.Read_GML('../0.preprocess_process/network_data/' + lang_code + 'wiki.gml')
    name_list = pd.read_pickle('../0.preprocess_process/network_data/'  + lang_code + 'wiki.pkl')
    target_nodes = pd.read_pickle('../2.language_link/target_files/target_set_' + lang_code + '.pkl')
    target_idx = name_list[name_list.isin(target_nodes)].index
    
    target_index = name_list[name_list == 'Science_and_technology_' + lang_code].index[0]
    G.vs['level'] = G.shortest_paths_dijkstra(source = target_index, mode= 'IN')[0]
    
    level_dict = {node.index:node['level'] for node in G.vs}
    adj_list = [G.neighbors(index, mode='OUT') for index in range(len(name_list))]
    
    total_length = len(adj_list)
    W_matrix = []
    for row in adj_list:
        dummy_dict = {}
        length = len(row)
        if length > 0:
            for node in row:
                dummy_dict[node] = np.power(k_factor, level_dict[node])
            dummy_dict = normalize(dummy_dict)     
        dummy_csr = scipy.sparse.csr_matrix((list(dummy_dict.values()), ([0]*len(dummy_dict),list(dummy_dict.keys()))), shape=(1, total_length))
        W_matrix.append(dummy_csr)
    W_matrix = scipy.sparse.vstack(W_matrix)
    
    splited_target_idx = np.array_split(target_idx, number_of_process)
    outputs = [Queue() for i in range(number_of_process)]  
    procs = [Process(target=calculate_ppr, args = (i, splited_target_idx[i], W_matrix, alpha, outputs[i])) for i in range(number_of_process)]
    
    for p in procs:
        p.start()
    print('result_collecting')
    
    result = {}
    for output in outputs:
        result.update(output.get())
    print('output_queue_closed')
    
    for output in outputs:
        output.close()
        
    print('process_closed')
    for p in procs:
        p.join()
    
    pickle.dump(result, open('ppr_result/' + lang_code + 'wiki_ppr.pkl+', 'wb'))