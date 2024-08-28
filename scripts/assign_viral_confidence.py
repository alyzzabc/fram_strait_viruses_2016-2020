import os
import pandas as pd
import numpy as np 

d = "/path/to/vs2-checkv/tables"

#Make list of combined vs2-checkv tables
with open("/path/to/list/of/filenames/" , "r") as tsv_list:
    for line in tsv_list: #m,n,p,q only for counting, but can be ignored
		m=0
		n=0
    	p=0
		q=0
    	tsv_file = os.path.normpath(d + line.strip()) #Create path to tsv file
    	df = pd.read_csv(tsv_file, sep='\t') #Read file and convert to dataframe
    	df = df.assign(screening_category='') #Add screening_category column to end of dataframe
    	#print(df)

    	#Read each row in tsv_file and assign category
    	for row in df.itertuples(): 
    		if row[15] > 0:
    			#print("keep1:" + row[1])
    			m=m+1
    			df.at[row.Index, "screening_category"] = "keep1"
    		elif (row[15] == 0 and row[16] == 0) or (row[15] == 0 and row[7] >= 0.95) or (row[15] == 0 and row[10] > 2):
    			#print("keep2:" + row[1])
    			n=n+1	
    			df.at[row.Index, "screening_category"] = "keep2"
    		elif (row[15] == 0 and row[16] == 1 and row[9] >=10000):
    			#print("manual_check:" + row[1])
    			p=p+1
    			df.at[row.Index, "screening_category"] = "manual_check"
    		else:
    			#print("discard:" + row[1])
    			q=q+1
    			df.at[row.Index, "screening_category"] = "discard"

    	print("Total number of contigs:" + str(m+n+p+q))
    	print("keep1:" + str(m))
    	print("keep2:" + str(n))
    	print("manual_check:" + str(p))
    	print("discard:" + str(q))
    	df.to_csv(d + line.strip() + "_screened.tsv", sep="\t", header=True, index=False, na_rep="NA") #Save assignments as tsv file
