import sys
from sys import stderr


"""
A load flow program
"""

def main():    
    str = "Calculate load flow for a network..."
    print(str)
    
    """
    Read in data for load flow studies:
    Data can be input from 
      local file system: csv file or excel file
      local database: MySQL, Oracle XE, or PostgreSQL
      Currently, implemented reading from a CSV file.
    
    Data will be processed to form B1P, B11P, etc 

    """     
    pqlfin()
   

    """
    Solve load flow        
    """
#    ///fdlf()


    """
    Output loadflow calculations to a file.
    """
#    ///pqlfout()



if __name__ == '__main__':
    main()
    

