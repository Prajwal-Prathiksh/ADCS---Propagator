tle_file = open('ins_2018.txt', 'r') 
# Reading The .txt file

tle = tle_file.readlines()

tle_2018 = []
# Storing the TLEs in a list, after suitable formatting.
for i in tle:
    tle_2018.append( i.rstrip() )   