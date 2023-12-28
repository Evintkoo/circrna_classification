from typing import List
import numpy as np
import math
import tqdm
import multiprocessing

NUCLEOTIDES_LIST = ["A", "G", "C", "U"]


def make_rna_combination_list(strand_length_target:int, cur_strand:str="", strand_cumulator:List[str]=[]) -> list:
    """Create list of combination of all possible RNA in the string

    Args:
        strand_length_target (int): _description_
        cur_strand (str, optional): _description_. Defaults to "".
        strand_cumulator (list, optional): _description_. Defaults to [].

    Returns:
        list: _description_
    """
    if cur_strand == "":
        strand_cumulator = []
    if len(cur_strand) == strand_length_target and cur_strand not in strand_cumulator:
        strand_cumulator.append(cur_strand)
    elif len(cur_strand) < strand_length_target:
        for nucleotides in NUCLEOTIDES_LIST:
            make_rna_combination_list(strand_length_target=strand_length_target, 
                                      cur_strand=cur_strand + nucleotides, 
                                      strand_cumulator=strand_cumulator)
    return strand_cumulator

def kmer_count(rna_strand:str, k_val:int, comb_list:List[str]) -> list:
    rna_strand = rna_strand.upper()

    nucleotides_count_cumulator = {}
    for i in range(k_val, len(rna_strand)):
        for nucleotides in comb_list:
            nucleotides_count_cumulator[nucleotides] = rna_strand.count(nucleotides)
    return nucleotides_count_cumulator

def rna_validator(rna_strand:str):
    rna_strand = rna_strand.upper()
    if set(rna_strand) - set("".join(NUCLEOTIDES_LIST)):
        return -1
    else:
        return 0 

def rna_to_dictionary(rna_strand: str, k_val:int, comb_list:List[str]):
    rna_strand = rna_strand.replace("t", "u")
    rna_strand = rna_strand.replace("T", "U")
    if rna_validator(rna_strand) == -1:
        print("invalid")
        return 0
    else:
        return kmer_count(rna_strand=rna_strand, k_val=k_val, comb_list=comb_list)
    
def calculate_k_value(rna_strand_list:List[str]):
    avg_length = np.mean([len(i) for i in rna_strand_list])
    return int(math.log(avg_length,len(NUCLEOTIDES_LIST)))
    
def rna_strand_list_to_matrix(rna_strand_list:List[str], k_val:int=None):
    if k_val == None:
        k_val = calculate_k_value(rna_strand_list)
    nucleotides_combination_list = make_rna_combination_list(strand_length_target = k_val)
    list_of_rna_dictionary = [rna_to_dictionary(rna_strand=rna_strand, k_val=k_val, comb_list=nucleotides_combination_list) for rna_strand in rna_strand_list if rna_to_dictionary(rna_strand=rna_strand, k_val=k_val, comb_list=nucleotides_combination_list) != 0]
    
    #for rna_strand in rna_strand_list:
    #    rna_dictionary = rna_to_dictionary(rna_strand=rna_strand, k_val=k_val, comb_list=nucleotides_combination_list)
    #    if rna_dictionary != 0:
    #        list_of_rna_dictionary.append(rna_dictionary)
    
    return list_of_rna_dictionary

class sequence2matrix:
    def __init__(self, k_val:int=None) -> None:
        self.k_value = k_val
    
    def map_func(self, rna_strand):
        return rna_to_dictionary(rna_strand=rna_strand, k_val=self.k_value, comb_list=self.rna_combination)
    
    def fit(self, list_of_rna_sequence:List[str], save_dir:str):
        if self.k_value == None:
            self.k_value = calculate_k_value(list_of_rna_sequence)
        self.rna_combination = make_rna_combination_list(strand_length_target=self.k_value)
        matrix = np.array([])
        num_workers = int(multiprocessing.cpu_count())
        chunksize = max(1,int(len(list_of_rna_sequence)/num_workers))
        
        pool = multiprocessing.Pool(num_workers)
        matrix = pool.map(self.map_func, list_of_rna_sequence)
        #for rna_strand in tqdm(list_of_rna_sequence,):
        #    vector = rna_to_dictionary(rna_strand=rna_strand, k_val=self.k_value, comb_list=self.rna_combination)
        #    if vector != 0:
        #        matrix.append(vector)
        #        np.save(save_dir, matrix)
        #np.save(save_dir, matrix)
        #np.array([rna_to_dictionary(rna_strand=rna_strand, k_val=self.k_value, comb_list=self.rna_combination) for rna_strand in list_of_rna_sequence if rna_to_dictionary(rna_strand=rna_strand, k_val=self.k_value, comb_list=self.rna_combination) != 0])
        np.save(save_dir+"_parameter", self.rna_combination)
        return np.array(matrix)