import numpy as np


#function that returns a dictionary having sequence name - sequence pairs
#download a multi-FASTA file from the NCBI database (optional: save it as seqdump.txt)

def cleaning(text):
    my_dict = {}
    with open(text, 'r') as file:
        for line in file:
            if line[0] == '>':
                try:
                    next_line = next(file)
                    my_dict[line.strip()] = next_line.strip()
                    temp = line.strip()
                except StopIteration:
                    print("end of the file")
            else:
                my_dict[temp] += line.strip()
    return my_dict

#function for counting a frequencies dictionary

def kf(seq, k):
    my_dic = {}
    for i in range(len(seq) - k+1):
            a = seq[i:i+k]
            if a not in my_dic:
                my_dic[a] = 1
            else:
                my_dic[a] += 1
    return my_dic

# dic = cleaning('seqdump.txt')
# for i, j in dic.items():
#     print(f'the freq for the sequence ({i}) is:\n')
#     print(f'{kf(j,3)}\n')
#outputs the k-mer frequencies for every sequence in the file 'seqdump.txt'


#function that calculates the cosine similarity of two given sequences in a dictionary format

def cossim(dic1, dic2):
    a = np.array(list(dic1.keys()))
    b = np.array(list(dic2.keys()))
    magnitude = np.linalg.norm(np.array(list(dic1.values())))*np.linalg.norm(np.array(list(dic2.values())))
    c = np.intersect1d(a,b)
    dot_pro = 0
    for m in c:
        dot_pro += dic1[m]*dic2[m]
    return (dot_pro/magnitude)

