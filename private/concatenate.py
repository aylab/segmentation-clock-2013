import os

def ensureDir(directory):
	if directory[-1] == '/':
		directory = directory[-1]
	if not os.path.exists(directory):
		return False
	return True

def main():
    print "HI"
    b = open("behavior_strict2.csv", "w")
    a = open("allpassed_strict2.csv", "w")
    folders = ["o0", "o3", "o5", "o6", "o7", "o8", "o11", "o12"]
        
    header = False
    
    parindex = 0
    # for each of the o folders
    #for folder in folders:
        # for each of the ox-i folders
        #for i in range(50):
        #       curdir = folder + "/" + folder + "-" + str(i)
        #       if ensureDir(curdir):
                    #for every line in the behavior file
    curb = open("behavior.csv")
    apass = open("allpassed.csv")
    
    hcurb = curb.readline()
    hpass = apass.readline()
    hpass = hpass.split(",")
    hpass = hpass[1:]
    hpass = ",".join(hpass)
    if not header:
        hstr = str(hcurb)
        hstr = hstr[:-1]
        hstr = hstr.split(",")
        hstr[4] = "per delta/wt"
        hstr[5] = "sync delta/wt"
        hstr[6] = "amp delta/wt"
        
        hstr[7] = "per her1/wt"
        hstr[8] = "sync her1/wt"
        hstr[9] = "amp her1/wt"
        
        hstr[10] = "per her7/wt"
        hstr[11] = "sync her7/wt"
        hstr[12] = "amp her7/wt"
        
        hstr[13] = "per her13/wt"
        hstr[14] = "sync her13/wt"
        hstr[15] = "amp her13/wt"
        
        hstr[16] = "per her713/wt"
        hstr[17] = "sync her713/wt"
        hstr[18] = "amp her713/wt\n"
        hcurb = ",".join(hstr)
        
        a.write(hpass)
        b.write(hcurb)
        header = True
    for line in curb:
        bline = line.split(",")
        bline[-1] = bline[-1][:-1]

        #if it shows passing behavior
        if len(bline) == 20:
            par = apass.readline()
            par = par.split(",")
            
            #parfolder = curdir + "/par" + par[0]
            #destfolder = "good sets/par" + str(parindex)
            #os.system('cp -r ' + parfolder + " " + destfolder)
            
            actualline = []
            actualline.extend(bline)
            actualline[4] = str(float(bline[5]) / float(bline[2])) #per delta/wt
            actualline[5] = str(float(bline[4]) / float(bline[1])) #sync delta/wt
            actualline[6] = str(float(bline[6]) / float(bline[3])) #amp delta/wt
            
            actualline[7] = str(float(bline[9]) / float(bline[2])) #per her1/wt
            actualline[8] = str(float(bline[15]) / float(bline[1])) #sync her1/wt
            actualline[9] = str(float(bline[10]) / float(bline[3])) #amp her1/wt
           
            actualline[10] = str(float(bline[11]) / float(bline[2])) #per her7/wt
            actualline[11] = str(float(bline[16]) / float(bline[1])) #sync her7/wt
            actualline[12] = str(float(bline[12]) / float(bline[3])) #amp her7/wt
        
            actualline[13] = str(float(bline[7]) / float(bline[2])) #per her13/wt
            actualline[14] = str(float(bline[17]) / float(bline[1])) #sync her13/wt
            actualline[15] = str(float(bline[8]) / float(bline[3])) #amp her13/wt
            
            actualline[16] = str(float(bline[13]) / float(bline[2])) #per her713/wt
            actualline[17] = str(float(bline[18]) / float(bline[1])) #sync her713/wt
            actualline[18] = str(float(bline[14]) / float(bline[3])) + "\n" #amp her713/wt
            actualline = actualline[:-1]
            if ((float(bline[16]) < 0.8 or float(actualline[12]) < 0.85) and float(actualline[1]) > 0.8 and float(bline[4]) < 0.7 and float(actualline[4]) > 1.04 and float(actualline[6]) < 0.9 and float(actualline[9]) < 1.06 and float(actualline[12]) < 1 and float(bline[17]) > 0.8 and float(actualline[15]) > 0.85 and float(bline[18]) > 0.8 and (float(bline[14]) / float(bline[3])) > 0.85):
                par[0] = str(parindex)
                par = par[1:]
                par = ",".join(par)
                actualline[0] = str(parindex)
                actualline = ",".join(actualline)
                b.write(actualline)
                a.write(par)
            
                parindex += 1
           
            

    curb.close()
    apass.close()

    a.close()
    b.close()

main()
