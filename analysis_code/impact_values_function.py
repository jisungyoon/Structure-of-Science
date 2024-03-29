from multiprocessing import Process, Queue
import numpy as np
import igraph as ig
import pandas as pd


def calculate_impact_values(G, lamda=1/2):
    in_list = G.get_adjlist(mode='IN')
    degree_in = np.array(G.degree(mode='IN'))
    out_list = G.get_adjlist(mode='OUT')
    degree_out = np.array(G.degree(mode='OUT'))
    number_of_node = len(degree_out)
    number_of_process = 30
    xs = np.array_split(np.arange(number_of_node), number_of_process)
    outputs = [Queue() for i in range(number_of_process)]   
    procs = [Process(target=calculate_edge_list, args = (i, in_list, out_list, degree_in, degree_out, xs[i], outputs[i])) for i in range(number_of_process)]
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

def calculate_edge_list(core_number,in_list, out_list, degree_in, degree_out, edge_range, output):
    print('CORE ' + str(core_number) + ' started')
    candidates_list = []
    lamda = 1/2
    epsilon = 0.0000001
    for index, i in enumerate(edge_range):
        candidates = {}
        if len(out_list[i]) == 1 :
            candidates[out_list[i][0]] = 1
        else:
            for j in out_list[i]:
                remain_node =  [x for x in in_list[j] if x is not int(i)]
                s_d_total = 0
                s_a_total = 0
                for k in remain_node:
                    if degree_in[k] != 0:
                        descendant_target_list = list(set(in_list[i]) & set(in_list[k]))
                        descendant_effet_values = [1/degree_out[l] for l in descendant_target_list]
                        s_d_total += np.nansum(descendant_effet_values)/degree_in[k]
                    if degree_out[k] != 0:
                        ancestor_target_list = list(set(out_list[i]) & set(out_list[k]))
                        ancestor_effect_values = [1/degree_in[l] for l in ancestor_target_list]
                        s_a_total += np.nansum(ancestor_effect_values)/degree_out[k]
                value = (s_d_total * lamda)  +  (s_a_total * (1-lamda))
                if value == 0:
                    value = epsilon
                candidates[j] = value
        candidates_list.append(candidates)
        if index % 1000 == 0:
            print("CORE " + str(core_number) +  ": " + str(index/len(edge_range)*100) + "%")
    output.put(candidates_list)
    print('sub_process_completed')
    
