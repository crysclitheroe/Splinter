"""

This script is written as the first in a set of simple tools to generate hypothetical mutant DNA/RNA sequences from a single wildtype sequence. mutagen.py takes as its input a single wildtype sequence and generates all the single_mutation permutations of that sequence. Ideally it reports the quantity of permutations and dumps all the hypothetical sequences to one of several filetypes for analysis.

Edit* the problem has changed slightly, only certain regions of the sequence need to be mutated. By happy accident the code works perfectly if the user inputs their wildtype sequence in lower case with the mutable regions of interest in uppercase. Yay for case-specific recognition. Now adding other information, to the data, such as actual mutation and position in the string, as well as mutation type.


"""
#import the libraries
import pandas as pd
import numpy as np
import re

#x-----------------------------------------------------------------------------------------------------------------------x
# All the Functions
#x-----------------------------------------------------------------------------------------------------------------------x

def single_mutants(wt,key):

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

   pattern = r'{}'.format(pos)
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



#x----------------------------------------------------------------------------------------------------------------------x

def library_maker(df):
 #call the single mutants function and transpose the df, in order to count right (i.e., along the index). need to figure out why I have to do this. there must be a way to to tell pandas to do x for every row

 dft = df.T
 shape = len(df.columns) #to get column index, which is the mutation info, which will be used to create the keys
 mutant_dictionary = {}
 
 #put sequences in the library
 for mu in dft:
  mutant_dictionary.update({df.ix[mu,(shape-1)] : df.ix[mu,0:(shape-1)].str.cat()})
 return mutant_dictionary

#x----------------------------------------------------------------------------------------------------------------------x

def double_mutants(library):
 
#essentially what we're doing here is looping through the single mutant library, calling the single_mutant function again on each single_mutant to make the double mutants (this will naturally create duplicates). serially append each sub df to make double mutant df, then drop duplicates.
 ddf = pd.DataFrame()
 shape = len(ddf.columns)

 for key, seq in library.items():
  ddf = ddf.append(single_mutants(seq,key), ignore_index = True)
 
 #and rows where no mutation info was set, i.e where the same mutation has already happened
 ddf = ddf[ddf.Mutation != '']


 print(ddf)
  
 #this drop duplicates thing was a pain to figure out - remember this!
 return ddf.drop_duplicates(ddf.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #return ddf.reset_index(drop = True)
                
#x-----------------------------------------------------------------------------------------------------------------------x 
#Function to create files

def mk_filer(df, wt, filename):

 dft = df.T
 shape = len(df.columns)

 #the fasta file writer is broken. Which is probably a good thin because I need to modernify it
 #open and write fasta file
 #from sys import argv
 #fafile = open('{}.fa'.format(filename), 'w')
 
 #remove the data to sequence strings
 #for mu in dft:

  #fafile.write('>mutant_{}'.format(df.Mutation))
  #fafile.write('\n')
  #fafile.write(df.ix[mu,0:(shape-1)].str.cat().upper())
  #fafile.write('\n')

 #fafile.close()

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
#Function to claculate the Single deletion mutant list

def del1_mutants(wt, pos):

 #start dictionary
 del1_diction = {}
 for ind, char in enumerate(wt):

  pattern = r'{}'.format(ind)
  matchObj = re.search(pattern, wt)

  print(matchObj)

  if char.isupper() == True and matchObj == None:
   key = '{}{}, {}'.format(wt[ind], ind+1, pos)
   mut = list(wt)
   mut[ind] = "_"
   del1_diction[key] = "".join(mut)

 return del1_diction

def sdel_mutants(libr):
 my_list = []
 table = []

 #had to do some interesting things with extend and append to get the array right
 for k, v in libr.items():
  my_list = []
  my_list.extend(list(v))
  my_list.append(k)
  table.append(my_list)
 
 del1 = pd.DataFrame(table)

 #stack overflow approach to solving the renaming thing
 new_col = int(len(del1.columns))
 names = del1.columns.tolist()
 names[names.index(new_col-1)] = 'Mutation'
 del1.columns = names
 print(del1)
 shape = len(del1.columns)
 #return del1.drop_duplicates(del1.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 return del1.reset_index(drop = True)

#x-----------------------------------------------------------------------------------------------------------------------x
#Function to claculate the Double deletion mutant list

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

  #stack overflow approach to solving the renaming mutation column for now
  new_col = int(len(temp_df.columns))
  names = temp_df.columns.tolist()
  names[names.index(new_col-1)] = 'Mutation'
  temp_df.columns = names
  
  del2_df = del2_df.append(temp_df, ignore_index = True)
  print(del2_df)
  shape = len(del2_df.columns)
 
 return del2_df.drop_duplicates(del2_df.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 #return del2_df.reset_index(drop = True)
 
#x-----------------------------------------------------------------------------------------------------------------------x
#Function to claculate the Triple deletion mutant list

def del3_mutants(del2s):
 del3_df = del2_mutants(del2s)
 shape = len(del3_df.columns)
 print(del3_df)
 #return del3_df.drop_duplicates(del3_df.columns[0:(shape-1)], keep = 'last').reset_index(drop = True)
 return del3_df.reset_index(drop = True)

#x-----------------------------------------------------------------------------------------------------------------------x
# a Function to place all mutants created into a dataframe

def full_mutant_df(smd,dmd,del1d,del2,del3):
 table = []
 for k, v in del1d.items():
  table.append(list(v))
 del1 = pd.DataFrame(table)
 frames = [smd,dmd,del1,del2,del3]
  
 sum_df = pd.concat(frames, ignore_index=True)
 sum_df = sum_df.fillna('')
 sum_df = sum_df.drop_duplicates()
 
 return sum_df.reset_index(drop = True)
 
#x-----------------------------------------------------------------------------------------------------------------------x
# Main - menu
#x-----------------------------------------------------------------------------------------------------------------------x

#First generate a single sequence for pilot 
user = input("\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n1) Single Deletion Mutant Library\n2) Double Deletion Mutant Library\n3) Triple Deletion Mutant Library\n4) Single Replacement Mutant Library\n5) Double Replacement Mutant Library\n6) Generate All Mutant Libraries\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\nEnter value to generate desired mutant library:  ").strip(" ")

while user != '':
 if user == "1":
  wild = input("\n\nEnter sequence to be mutated here, with the mutable region in upper case: \n").strip(" ")
  #calculate permutation number here
  filename = input("Enter a file name for the library (no spaces or other punctuation or special characters - please!): ").strip(" ")
  print("\n\nA single deletion mutant library of xxx sequences is being generated...\n")
  mk_filer(single_mutants(wild, ''), wild, filename)
  user = input("\n\n...Complete!\nFile xxx has been created in the working directory.\n\nPress Enter to return to Main Menu or Escape to exit.\n\n")
  # a is to generate a list of single deletion mutants, need to show permutation calculation number and give a summary of the created file
 elif user == "2":
  break # b is to generate a list of double deletion mutants, need to show permutation calculation number and give a summary of the created file
 elif user == "3":
  break # 3 is to generate a list of Triple deletion mutants, need to show permutation calculation number and give a summary of the created file
 elif user == "4":
  break # 4 is to generate a list of single replacement mutants, need to show permutation calculation number and give a summary of the created file
 elif user == "5":
  break # 5 is to generate a list of double replacement mutants, need to show permutation calculation number and give a summary of the created file
 elif user == "6":
  break # 6 is to generate a list of all replacement and deletion mutants, need to show number and give a summary of the created file
 else:
  user = input("\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n1) Single Deletion Mutant Library\n2) Double Deletion Mutant Library\n3) Triple Deletion Mutant Library\n4) Single Replacement Mutant Library\n5) Double Replacement Mutant Library\n6) Generate All Mutant Libraries\n\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx WELCOME TO SPLINTER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n\nEnter value to generate desired mutant library:  ").strip(" ")

 
#user = input("Enter sequence to be mutated here, with the mutable region in upper case: ").strip(" ")
#user = "gagaccgcaactgaaaagttgtTATCACTTGTCGTAAGACACTTTGGATGggttGAAgttctgcacgtagaagcaaaaggc" 
#user = 'agAGCTct'


#
#mk_filer(double_mutants(library_maker(single_mutants(user,''))),user, "Double_Mutant_Library")
#mk_filer(all_mu, user, "Mutant_Library")
#mk_filer(del2_mutants(del1_mutants(user,'')),user,'DblDel_Mutant_Library')
#mk_filer(del3_mutants(library_maker(del2_mutants(del1_mutants(user,'')))),user,'TrpDel_Mutant_Library')
#mk_filer(sdel_mutants(del1_mutants(user,'')),user,'SglDel_Mutant_Library')

#all_mu = full_mutant_df(single_mutants(user),double_mutants(library_maker(single_mutants(user))),del1_mutants(user),del2_mutants(del1_mutants(user)),del3_mutants(library_maker(del2_mutants(del1_mutants(user)))))





