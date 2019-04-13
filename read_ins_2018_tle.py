tle_file = open('ins_2018.txt', 'r') 
# Reading The .txt file
tle = tle_file.readlines()
tle_2018 = []
# Storing the TLEs in a list, after suitable formatting.
for i in tle:
    tle_2018.append( i.rstrip() ) 
tle_file.close()


tle_file = open('pratham_tle.txt', 'r')
# Reading The .txt file
tle = tle_file.readlines()
pratham_tle_2018 = []
# Storing the TLEs in a list, after suitable formatting.
for i in tle:
    pratham_tle_2018.append( i.rstrip() ) 
tle_file.close()

tle_file = open('pisat_tle.txt', 'r')
# Reading The .txt file
tle = tle_file.readlines()
pisat_tle_2018 = []
# Storing the TLEs in a list, after suitable formatting.
for i in tle:
    pisat_tle_2018.append( i.rstrip() ) 
tle_file.close()

tle_file = open('niusat_tle.txt', 'r')
# Reading The .txt file
tle = tle_file.readlines()
niusat_tle_2018 = []
# Storing the TLEs in a list, after suitable formatting.
for i in tle:
    niusat_tle_2018.append( i.rstrip() ) 
tle_file.close()
