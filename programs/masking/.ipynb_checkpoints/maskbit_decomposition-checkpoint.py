'''
BITMASK DICTIONARY:
https://www.legacysurvey.org/dr9/bitmasks/
'''

import os
homedir = os.getenv("HOME")

import numpy as np

#create dictionary with exponents (n) and their resulting 2**n value
def create_power_dictionary(min_power=0, max_power=12):

    two_exp = np.arange(min_power,max_power+1,1)
    two_array = two_exp.copy()
    
    for i in two_array:
        two_array[i] = 2**i
        
    powers_of_two = dict(zip(two_exp,two_array))
    
    return(powers_of_two)


#decompose every pixel of the r-band mask into a list of the constituent powers
#this function is run for EVERY pixel. yeesh. :-(
def power_decomposition(old_pixel_val, power_dict):
    
    #create empty decomposition list (which will contain the exponents that comprise the mask)
    decomposition = []
    
    #isolate old pixel value
    n = old_pixel_val
    
    #loop through every exponent in power dictionary
    for i in range(len(power_dict)):
        
        #begin with largest power of 2
        dict_index = len(powers_of_two)-(1+i)
        
        #if n is already zero, finish
        if n == 0:
            break
        
        #if n is *exactly* equal to a power of 2, append to decomposition list and finish
        if n == powers_of_two[dict_index]:
            decomposition.append(dict_index)
            break
        
        #begin with largest power of 2 remaining
        dict_index = len(powers_of_two)-(1+i)
        #determine the number to subtract from n
        byte_to_subtract = powers_of_two[dict_index]

        #in the case where the number to subtract is larger than n, skip to the next power of 2
        if (n-byte_to_subtract)<0:
            pass
        
        #else, subtract number from n and append the corresponding power of 2 to the decomposition list
        if (n-byte_to_subtract)>=0:
            n -= byte_to_subtract
            decomposition.append(dict_index)
    
    #return a list of every powers-of-2 exponent contributing to the maskbit pixel value
    return(decomposition)