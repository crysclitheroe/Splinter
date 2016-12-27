"""

This script is written as the first in a set of simple tools to generate hypothetical mutant DNA/RNA sequences from a single wildtype sequence. mutagen.py takes as its input a single wildtype sequence and generates all the single_mutation permutations of that sequence. Ideally it reports the quantity of permutations and dumps all the hypothetical sequences to excel, as this is commonly the format that one sends the library to sequence generation companies (e.g. Applied Biosciences).


"""
#import the libraries
import pandas as pd
import numpy as np
import re

#x-----------------------------------------------------------------------------------------------------------------------x
# All the Functions
#x-----------------------------------------------------------------------------------------------------------------------x

#x-----------------------------------------------------------------------------------------------------------------------x
# Option 1) Function to calculate the Single deletion mutant list 

def del1_mutants(wt, pos):

 #start dictionary
 del1_diction = {}
 for ind, char in enumerate(wt):

  #start the pattern for recognising and eluding positions that have been deleted 
  pattern = r'{}'.format(ind)
  matchObj = re.search(pattern, wt)

  #if index is in the mutable region and the mutation is not repeated
  if char.isupper() == True and matchObj == None:
   key = '{}{}, {}'.format(wt[ind], ind+1, pos)
   mut = list(wt)
   #Prof's feature request, place a marker whe deletion has occurred
   mut[ind] = "_"
   del1_diction[key] = "".join(mut)

 #rturn a library that later functions can work with
 return del1_diction

def sdel_mutants(libr):
 my_list = []
 table = []

 #rearrange the library to get the array structure consistent
 for k, v in libr.items():
  my_list = []
  my_list.extend(list(v))
  my_list.append(k)
  table.append(my_list)
 #send to a df
 del1 = pd.DataFrame(table)

 #Rename the column with Mutation info
 new_col = int(len(del1.columns))
 names = del1.columns.tolist()
 names[names.index(new_col-1)] = 'Mutation'
 del1.columns = names
 
 #get df attributes to process and return only unique sequence data
 shape = len(del1.columns)
 return del1.drop_duplicates(del1.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #Comment out above and uncomment below if you are curious about how your data looks with duplicates!
 #return del1.reset_index(drop = True)

#x-----------------------------------------------------------------------------------------------------------------------x
# Option 2) Function to calculate the Double deletion mutant list

def del2_mutants(del1s):

 del2_df = pd.DataFrame()

 for key, mut in del1s.items():
  temp = del1_mutants(mut,key) 
  table = []

  for k, v in temp.items():

   pattern = r'{}'.format(key)
   matchObj = re.search(pattern, mut)

   my_list = []
   my_list.extend(list(v))
   my_list.append(k)
   table.append(my_list)
  temp_df = pd.DataFrame(table)

  #Renaming the mutation column 
  new_col = int(len(temp_df.columns))
  names = temp_df.columns.tolist()
  names[names.index(new_col-1)] = 'Mutation'
  temp_df.columns = names
  
  #get df attributes to process and return only unique sequence data
  del2_df = del2_df.append(temp_df, ignore_index = True)
  shape = len(del2_df.columns)
 
 return del2_df.drop_duplicates(del2_df.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #Comment out above and uncomment below if you are curious about how your data looks with duplicates!
 #return del2_df.reset_index(drop = True)
 
#x-----------------------------------------------------------------------------------------------------------------------x
# Option 3) Function to calculate the Triple deletion mutant list

def del3_mutants(del2s):
 del3_df = del2_mutants(del2s)
 shape = len(del3_df.columns)
 
 return del3_df.drop_duplicates(del3_df.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #Comment out above and uncomment below if you are curious about how your data looks with duplicates!
 #return del3_df.reset_index(drop = True) 

#x-----------------------------------------------------------------------------------------------------------------------x
#Functions to calculate number of deletion mutant combinations, i.e. size of the library, this is a kind of internal check

def factor(val):

 factorial = 1
 for i in range(1,val+1):
  factorial = factorial*i
 return factorial

def del_combinator(wild, dels):
 
 #based on this calculator: http://www.calculatorsoup.com/calculators/discretemathematics/combinations.php
 mr = sum([c.isupper() for c in wild])
 combos = factor(mr) / (factor(dels)*factor(mr-dels))
 return combos

#x-----------------------------------------------------------------------------------------------------------------------x
# Option 4) Functions to calculate a single replacement mutant library

def single_mutants(wt, key):

 #Read in the sequence and send to a dataframe
 table = []
 table.append(list(wt))
 df = pd.DataFrame(table)
 mut_colno = len(wt) 

 #prepare a list of available bases
 bases = ['A', 'C', 'G', 'T']
 altbase = []

 #set the mutant seq number counter
 seq = 0
 counter = 0

 #iterate through the base list, until we hit a match with bases of the wildtype sequence (df columns)
 for base in bases:
  for pos in df.columns:

   #create regext of current position, to search mutation info, preventing back-mutations)
   pattern = r'{}'.format(pos+1)
   matchObj = re.search(pattern, key)

   #so if the base matches the wildtype and the position of previous mutation has not been covered
   if base == df.ix[0,pos] and matchObj == None:
    
    #prepare a list of alternate bases and use those to start new corresponding mutant seqs
    altbase.extend(bases)
    altbase.remove(base)
    df.set_value(seq+1, pos, altbase.pop())
    df.set_value(seq+2, pos, altbase.pop())
    df.set_value(seq+3, pos, altbase.pop())

    #summarise mutation info to new column, provided the mutants have not been seen before
    df.set_value(seq+1, mut_colno, '{}{}{}, {}'.format(base,pos+1,df.ix[seq+1,pos],key))
    df.set_value(seq+2, mut_colno, '{}{}{}, {}'.format(base,pos+1,df.ix[seq+2,pos],key))
    df.set_value(seq+3, mut_colno, '{}{}{}, {}'.format(base,pos+1,df.ix[seq+3,pos],key))
    seq += 3
 
 #replace the rest of the dataframe (empty mutant seq bases) with wildtype bases
 df = df.fillna(df.ix[0,df.columns])
 
 #change the column name for the mutations
 df = df.rename(columns = {mut_colno: 'Mutation'})

 #drop the wild type and rows where no mutation info was set, i.e where the same mutation has already happened
 df = df[df.Mutation.notnull()]

 return df.reset_index(drop = True)

#x----------------------------------------------------------------------------------------------------------------------x
# An in-between step that I should probably see if I can remove at some point

def library_maker(df):
 #call the single mutants function and transpose the df, in order to count right (i.e., along the index)
 dft = df.T
 shape = len(df.columns) #to get column index, which is the mutation info, which will be used to create the keys
 mutant_dictionary = {}
 
 #put sequences in the library for the double replacement mutant library
 for mu in dft:
  mutant_dictionary.update({df.ix[mu,(shape-1)] : df.ix[mu,0:(shape-1)].str.cat()})
 
 return mutant_dictionary

#x----------------------------------------------------------------------------------------------------------------------x
# Option 5) Functions to calculate a double replacement mutant library

def double_mutants(library):
 
#essentially what we're doing here is looping through the single mutant library, calling the single_mutant function again on each single_mutant to make the double mutants (this will naturally create duplicates). serially append each sub df to make double mutant df, then drop duplicates.
 ddf = pd.DataFrame()
 shape = len(ddf.columns)

 for key, seq in library.items():
  ddf = ddf.append(single_mutants(seq, key), ignore_index = True)
 
 return ddf.drop_duplicates(ddf.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #Comment out above and uncomment below if you are curious about how your data looks with duplicates!
 #return ddf.reset_index(drop = True)
                
#x-----------------------------------------------------------------------------------------------------------------------x 
# Function to create excel files

def mk_filer(df, wt, filename):

 dft = df.T
 shape = len(df.columns)

 #open and write to excel
 import xlwt
 book = xlwt.Workbook(encoding='utf-8')
 sh = book.add_sheet(filename)

 #create headers
 sh.write(0,0,'ID')
 sh.write(0,1,'Sequence')
 sh.write(1,0,'Wildtype')
 sh.write(0,2,'Mutation')
 sh.write(1,1,wt)

 for mu in dft:
  r=int(mu)+2
  sh.write(r,0,int(mu))
  sh.write(r,1,df.ix[mu,0:(shape-1)].str.cat().upper())
  sh.write(r,2,df.Mutation[mu])

 book.save('{}.xls'.format(filename))

#x-----------------------------------------------------------------------------------------------------------------------x
# Main - menu
#x-----------------------------------------------------------------------------------------------------------------------x

#First generate a single sequence for pilot 
user = input("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n1) Single Deletion Mutant Library\n2) Double Deletion Mutant Library\n3) Triple Deletion Mutant Library\n4) Single Replacement Mutant Library\n5) Double Replacement Mutant Library\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\nEnter value to generate desired mutant library:  ").strip(" ")
wild = ''

while user != '':
 if user == "1":
  if wild == '':
   wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
  filename = input("Enter a file name for the library (no spaces, punctuation or special characters - please!): ").strip(" ")
  print("\n\nA single deletion mutant library of {} unique sequences is being generated...\n".format(del_combinator(wild, 1)))
  singledf = sdel_mutants(del1_mutants(wild,''))
  mk_filer(singledf, wild, filename)
  user = input("...Complete!\nFile {}.xls has been created in the working directory.\n\nPress r, then Enter to return to Main Menu or Enter to exit.\n\n".format(filename))

 elif user == "2":
  if wild == '':
   wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
  filename = input("Enter a file name for the library (no spaces, punctuation or special characters - please!): ").strip(" ")
  print("\n\nA double deletion mutant library of {} unique sequences is being generated...\n".format(del_combinator(wild, 2)))
  doubledf = del2_mutants(del1_mutants(wild, ''))
  mk_filer(doubledf, wild, filename)
  user = input("...Complete!\nFile {}.xls has been created in the working directory.\n\nPress r, then Enter to return to Main Menu or Enter to exit.\n\n".format(filename))

 elif user == "3":
  if wild == '':
   wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
  filename = input("Enter a file name for the library (no spaces, punctuation or special characters - please!): ").strip(" ")
  print("\n\nA triple deletion mutant library of {} unique sequences is being generated...\n".format(del_combinator(wild, 3)))
  tripledf = del3_mutants(library_maker(del2_mutants(del1_mutants(wild,''))))
  mk_filer(tripledf, wild, filename)
  user = input("...Complete!\nFile {}.xls has been created in the working directory.\n\nPress r, then Enter to return to Main Menu or Enter to exit.\n\n".format(filename))

 elif user == "4":
  if wild == '':
   wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
  filename = input("Enter a file name for the library (no spaces, punctuation or special characters - please!): ").strip(" ")
  srepdf = single_mutants(wild,'')
  mk_filer(srepdf, wild, filename)
  user = input("...Complete!\nFile {}.xls has been created in the working directory.\n\nPress r, then Enter to return to Main Menu or Enter to exit.\n\n".format(filename))

 elif user == "5":
  if wild == '':
   wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
  filename = input("Enter a file name for the library (no spaces, punctuation or special characters - please!): ").strip(" ")
  drepdf = double_mutants(library_maker(single_mutants(wild,'')))
  mk_filer(drepdf, wild, filename)
  user = input("...Complete!\nFile {}.xls has been created in the working directory.\n\nPress r, then Enter to return to Main Menu or Enter to exit.\n\n".format(filename))

 else:
  user = input("\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n1) Single Deletion Mutant Library\n2) Double Deletion Mutant Library\n3) Triple Deletion Mutant Library\n4) Single Replacement Mutant Library\n5) Double Replacement Mutant Library\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\nEnter value to generate desired mutant library:  ").strip(" ")

 






