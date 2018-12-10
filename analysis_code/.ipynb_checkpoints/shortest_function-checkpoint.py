from multiprocessing import Process, Queue
import numpy as np
import igraph as ig
import pandas as pd

def calculate_all_shortest_path(G, target_index):
    print('initiating calculation')
    number_of_node = len(G.vs)
    number_of_process = 20
    xs = np.array_split(np.arange(number_of_node), number_of_process)
    outputs = [Queue() for i in range(number_of_process)]   
    procs = [Process(target=calculate_shortest_path, args = (i, G, xs[i], target_index, outputs[i])) for i in range(number_of_process)]
    for p in procs:
        p.start()
    print('result_collecting')
    result = []
    result_candidates = []
    for output in outputs:
        result += output.get()
    print('output_queue_closed')
    for output in outputs:
        output.close()
    print('process_closed')
    for p in procs:
        p.join()
        
    return result
        

def calculate_shortest_path(core_number, graph, node_range, target_index, output):
    result_array = []
    for index, i in enumerate(node_range):
        result_array.append(graph.get_all_shortest_paths(i, target_index, mode='OUT')) 
        if index % 1000 == 0:
            print("CORE " + str(core_number) +  ": " + str(index/len(node_range)*100) + "%")
    output.put(result_array)        
    print('core ' + str(core_number) + ' process_completed')   