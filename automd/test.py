import os

if __name__ == "__main__" :

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
