filenames = ['24_10_07/zpts.log','24_10_14/zpts.log','24_10_17/zpts.log','24_10_19/zpts.log','24_10_21/zpts.log']
with open('zeros.txt','w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())
