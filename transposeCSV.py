#!/usr/bin/python
#
# 2014-03-03 Shane Caldwell
#   Imports a .csv from Excel and transposes it
#   The input file is formatted for easy editing.
#   The output file is formatted for feeding into C structs.
#

import sys, csv

def main():
  filename = sys.argv[1]
  lis = []
  with open(filename,'r') as infile:
    with open(filename+'_transposed','w') as outfile:
      reader = csv.reader(infile, delimiter=',')#dialect='excel')
      writer = csv.writer(outfile, delimiter=',')
      for row in reader:
        #print len(row)
        #print row
        lis.append(row)
      #print lis
      #lis = zip(*lis)
      #print lis
      #for case in zip(*lis):
      #  writer.writerow(case)
      writer.writerows(zip(*lis))
      
#  infile  = open(sys.argv[1],'r')
#  outfile = file(sys.argv[1]+'_transposed','w')
#  #with open('bdn_cases.csv') as csvFile
#  lis=[row.split(',') for row in infile]
#  
#  #print lis
#  #print zip(*lis)
#  
#  for case in zip(*lis):
#    print case
#    print '*' * 80
#    for item in case:
      #print(entry+',')
      #outfile.write(item+',')
    #print('\n')
    #outfile.write('\n')
  
#  infile.close()
#  outfile.close()
  
if __name__ == '__main__':
  main()
