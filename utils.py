from typing import List

def cluster(l:List[float])->List[List[float]]: #cluster together contiguous numbers in a list [1,2,5,6] -> [[1,2],[5,6]]
    res = [[l[0]]]
    for i in l:
        if i == res[-1][-1] + 1:
            res[-1].append(i)
        else:
            res.append([i])
    return res[1:]

#split a list into parts dictated by the elements of the parts array
# split_lst([1,2,3,4,5,6,7,8], [3,5]) -> [[1, 2, 3], [4, 5, 6, 7, 8]]
def split_lst(lst:List[float],parts:List[int])->List[List[float]]: 
    res = []
    i,j = 0,0
    while j < len(parts):
        res.append(lst[i:i+parts[j]])
        i, j = i + parts[j], j + 1
    return res