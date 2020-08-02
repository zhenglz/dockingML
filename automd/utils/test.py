def passwd() :

    # cd pwd
    os.chdir(os.getcwd())

    bkp = "bkp.passwd"
    pwd = "passwd"

    loginList = []

    if os.path.exists(bkp) and os.path.exists(pwd) :

        bkp_lines = open(bkp).readlines()
        pwd_lines = open(pwd).readlines()

        if bkp_lines == pwd_lines :
            pass

        else :
            with open(pwd) as lines :

                tofile = open(bkp, 'w')
                for s in lines :
                    user = s.split(":")[0]

                    if "homes" in s and user not in loginList :
                        loginList.append(user)

                    if user in loginList :
                        infor = s.split(":")
                        infor[-1] = "/bin/sh"
                        tofile.write(":".join(infor) + " \n")
                    else :
                        tofile.write(s)
                tofile.close()
    else :
        os.system("touch %s"%bkp)

if __name__ == "__main__" :

    from automd.utils.gentop import *

    mol2File = sys.argv[1]
    fileList = []

    with open(mol2File) as lines :
        for s in [ x for x in lines if "#" not in x ] :
            fileList.append(s.split()[0])
    #start, end = int(sys.argv[2]), int(sys.argv[3])

    ## split mol2 file
    #mol = autoMD.ExtractPDB()
    #mol.extract_all(mol2File, "M")
    #nomol = len(mol.indexCoord(mol2File))

    #if end >= nomol :
    #    end = nomol
    AMBERHOME = "/home/liangzhen/software/amber14/"

    os.system("export AMBERHOME=/home/liangzhen/software/amber14/")

    for i in range(len(fileList)) :
        inpcrd = fileList[i]

        if not os.path.exists('./M'+str(i)) :
            os.mkdir('./M'+str(i))
        os.system("mv %s ./M%d/ " % (inpcrd, i))
        os.chdir('./M'+str(i))

        ## run antechamber
        os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
        top = GenerateTop()
        sh = top.runAntechamber(inpcrd, netCharge=False)
        job = sp.Popen("sh %s %s %s" % (sh, inpcrd, "M"+str(i)), shell=True)
        job.communicate()

        top.runObabel('obabel', inpcrd, "M"+str(i)+".pdb")
        structure = "M"+str(i)+".pdb"

        top.gmxTopBuilder("frcmod." + "M" + str(i),
                          "prep." + "M" + str(i),
                          structure, "M" + str(i),
                          None, 0, None, 0, FField=["gaff"])

        #top.removeMolInfor("M" + str(i))

        os.chdir("../../")
