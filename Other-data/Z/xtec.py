
#!/usr/bin/env python

import os

# loop over all folders & subfolders, identify all candidate files
for root, dirs, files in os.walk("./CT data"):

    for file in files:
        if file.endswith(".xtec"):
             print(os.path.join(root, file))

for root, dirs, files in os.walk("./CT data"):

    for file in files:
        if file.endswith(".xtekct"):
             print(os.path.join(root, file))



print("Bloops")

newfile = open("newfile.txt", "w")

# loop over all folders, read file contents, write to new file
for root, dirs, files in os.walk("./CT data"):
    for file in files:
        if file.endswith(".xtekct"):
		with open (os.path.join(root, file), "r") as myfile:
		    mytext = myfile.readlines()
		    mytext = ','.join(map(str, mytext)) 
		    newfile.write(mytext)
		    print(mytext)

newfile.close()


print("******************************************************")

#import csv

#with open("output.csv", "wb") as csv_file:
 #   writer = csv.writer(csv_file)
  #  for root, dirs, files in os.walk("./CT data"):
   # 	for file in files:
    #        if file.endswith(".xtekct"):
	#	with open (os.path.join(root, file), "r") as myfile:
	#	    mytext = myfile.readlines()
	#	    mytext = ','.join(map(str, mytext)) 
         #           writer.writerows(mytext)
