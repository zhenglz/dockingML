from .cmap import ContactMap
import os, sys
import urllib
from collections import OrderedDict
from os import environ
import subprocess as sp

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False

class FixPDB :
    def __init__(self):
        pass

    def pdbDownloader(self, pdbCode, pdbout):
        if not os.path.exists(pdbCode + '.pdb'):
            try :
                source = urllib.urlopen("http://www.rcsb.org/pdb/files/" + pdbCode + ".pdb")
                with open(pdbCode + '.pdb', 'wb') as target:
                    target.write(source.read())
            except urllib.URLError :
                print( "URL error")

        return(1)

    def saveHETATM(self, pdbin, chain=['A'], waterOnly=False):
        '''

        :param pdbin:
        :param chain:
        :param waterOnly: only save water molecule
        :return:
        '''
        with open(pdbin+"_HETATM", 'wb') as tofile :
            try :
                with open(pdbin) as lines :
                    for s in lines :
                        if len(s.split()) \
                                and s.split()[0] == "HETATM" \
                                and s[21] in chain :
                            if waterOnly and s[17:20].strip() in ["WAT","HOH"]:
                                tofile.write(s)
                            else :
                                tofile.write(s)
            except IOError :
                print("File %s not exist \n"% pdbin)

        return (pdbin + "_HETATM")

    def saveLigand(self, pdbin, chain, ligCode) :
        with open(pdbin + "_Ligand", 'wb') as tofile :
            try:
                with open(pdbin) as lines:
                    for s in lines:
                        if len(s.split()) \
                                and s.split()[0] in ["HETATM","ATOM"] \
                                and s[21] in chain\
                                and s[17:20].strip() == ligCode :
                            tofile.write(s)

            except IOError:
                print("File %s not exist \n" % pdbin)

        return (pdbin + "_Ligand")

    def checkCrash(self, pdbfiles, cutoff,
                   loopRange, ligandRange,
                   chainLoop, chainLigand,
                   maxCoordinate=10
                   ):
        '''
        check number of contacts between the newly modeled loop with the ligand
        if too many contacts, crashes are existed

        :param pdbfiles:
        :param cutoff:
        :param loopRange:
        :param ligandRange:
        :param chainLoop:
        :param chainLigand:
        :param maxCoordinate:
        :return: coordination number ( number of contacts)
        '''
        cmap = ContactMap(pdbfiles[0])
        coordinate = 0

        pdblist = []
        for item in pdbfiles :
            if '.pdb' not in item :
                if os.path.exists(item+'.pdb') :
                    pdblist.append(item+'.pdb')
                elif os.path.exists(item + '.pdbqt'):
                    pdblist.append(item + '.pdbqt')
                else :
                    print("File %s not found! "%item)
            else :
                pdblist.append(item)

        #try :
        coordinate = cmap.coordinationNumber(pdblist, cutoff,
                                            range(loopRange[0], loopRange[-1]+1),
                                            range(ligandRange[0], ligandRange[-1]+1),
                                            [chainLoop, chainLigand],
                                            ['heavy','all'],
                                            False, maxCoordinate
                                            )
        #except IOError:
            #print("Calculate coordination number between receptor and ligand failed!")

        return (coordinate)

    def addLoopLOOPY(self, pdbin,
                          loopRange,
                          loopSeq,
                          loopchain,
                          loopyexe='loopy',
                          num_mode = 10,
                          iteration=200,
                          ligandName='',
                          ligandNdx=1,
                          ligandchain='',
                          verbose=False
                          ):
        '''
        perform loop building with LOOPY
        LOOPY: https://honiglab.c2b2.columbia.edu/software/Jackal/Jackalmanual.htm#loopy
        Usage: ./loopy -o=900-908 -r=ALSVEFPEM -n=10 -i=200 2ovh.pdb -c=A > 2ovh_899_909

        Notes: in the process of loop modeling, the Ligand and the HETATM atom records are
        all lost. thus if you intend to keep them, you could save them ahead using the
        self.saveLigand and self.saveHETATM methods

        :param pdbin:
        :param loopRange:
        :param loopSeq:
        :param loopyexe:
        :param verbose:
        :return:
        '''

        cmd = "%s -o=%d-%d -r=%s -n=%d -i=%d -c=%s %s > %s_adding_loop_%d_%d.log" % \
              (
                  loopyexe, loopRange[0], loopRange[1], loopSeq,
                  num_mode, iteration,
                  loopchain, pdbin, pdbin,
                  loopRange[0], loopRange[1]
              )
        if not verbose :
            cmd += " > /dev/null 2>&1 "

        output = sp.check_output(cmd, shell=True)
        finalModel = pdbin + "_looper_0.pdb"
        if verbose :
            print(output)

        model_crash = OrderedDict()
        for i in range(num_mode) :
            model =  pdbin + "_looper_" + str(i) + ".pdb"
            # check crash. Search for any crash for the loop and the ligand.
            crash = self.checkCrash([model, pdbin], 1.5, loopRange, [ligandNdx,], loopchain, ligandchain, 10)
            model_crash[model] = crash

        #model_crash = OrderedDict(sorted(model_crash.items(), key=lambda x: x[1], reverse=False))
        #crashIndx = 0
        key = pdbin + "_looper_" + str(i) + ".pdb"
        for key in model_crash.keys():
            if model_crash[key] <= 1.0 :
                finalModel = key
                break
            else :
                pass
        if key == len(model_crash.keys()[-1]) :
            model_crash = OrderedDict(sorted(model_crash.items(), key=lambda x: x[1], reverse=False))
            finalModel = model_crash.keys()[0]

        if verbose :
            print("Checking crash in the LOOPY models ")
            print(model_crash)

        if model_crash.values()[0] > 10.0 :
            print("Warning: the LOOPY model %s is not perfect, crashes occur! ")

        #finalModel = model_crash.keys()[0]

        return (finalModel)

    def removeRegions(self, filename, residuesNdx, chain, pdbout="temp_1.pdb", verbose=False):
        '''
        input a pdbfile, remove the selected residues

        :param filename: input pdbfile
        :param residuesNdx: the range of residue index for deleting
        :param chain: which chain to perform delete
        :param pdbout: the output pdb file name
        :return: 1
        '''
        tofile = open(pdbout,'w')
        with open(filename) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) in residuesNdx :
                    pass
                else :
                    tofile.write(s)
        tofile.close()

        return(1)

    def addModeledRegions(self, basepbd, modeledpdb,
                          modelNdx, removeNdx, chain,
                          pdbout="temp.pdb",
                          verbose=False
                          ):
        tofile = open(pdbout, 'w')
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) < removeNdx[0] :
                    tofile.write(s)
                else :
                    pass

        ## addd modeled pdb here
        with open(modeledpdb) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and int(s[22:26].strip()) in modelNdx :
                    tofile.write(s)
                else :
                    pass

        ## add the following part of the original PDB
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and s[21] == chain and int(s[22:26].strip()) > removeNdx[-1]:
                    tofile.write(s)
                else :
                    pass

        tofile.close()
        return(1)

    def addhydrogenReduce(self, pdbin, pdbout='outputH.pdb',reduce='reduce', flipH=True,verbose=True):
        if flipH :
            cmd = "%s -FLIP %s > %s " % (reduce, pdbin, pdbout)
        else :
            cmd = "%s -NOFLIP %s > %s " % (reduce, pdbin, pdbout)
        job = sp.check_output(cmd, shell=True)

        if verbose :
            print(job)

        return 1

    def addMissingAtoms(self, pdbin, pdb2pqrPy,
                        ff='amber', pqrout='output.pqr',
                        addhydrogen=True,
                        ph=7.0, verbose=False ):
        '''
        add missing heavy atoms, and assign hydrogens according to the ph

        :param pdbin:
        :param pdb2pqrPy:
        :param ff:
        :param pqrout:
        :param addhydrogen:
        :param ph:
        :return:
        '''
        if addhydrogen :
            cmd = "python %s -ff %s --with-ph %f --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ph, ff+"_"+pqrout, pdbin, pqrout)
        else :
            cmd = "python %s -ff %s --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ff+"_"+pqrout, pdbin, pqrout)
        job = sp.Popen(cmd, shell=True)
        job.communicate()
        job.kill()

        return(ff+"_"+pqrout)

    def addLoopsAutoModel(self, pdbCode, chain1,
                         alignCode, chain2,
                         loopRange,
                         verbose=False):
        '''
        model missing loop region
        :param pdbCode:
        :param chain1:
        :param alignCode:
        :param chain2:
        :param loopRange:
        :param verbose: print more information
        :return: the final model pdb file
        '''
        if MODELLER_EXIST :
            if verbose :
                print("MODELLER exist, perform alignment and loop refinement.")
                log.verbose()

            for pdb in [pdbCode]+[alignCode] :
                if not os.path.exists(pdb+'.pdb') :
                    self.pdbDownloader(pdb, pdb+'.pdb')

            env = environ()
            # directories for input atom files
            env.io.atom_files_directory = ['.', '../atom_files']

            aln = alignment(env)
            mdl_1 = model(env, file=pdbCode, model_segment=('FIRST:%s'%chain1, 'LAST:%s'%chain1))
            aln.append_model(mdl_1, align_codes=pdbCode, atom_files=pdbCode+'.pdb')
            mdl_2 = model(env, file=alignCode, model_segment=('FIRST:%s'%chain2, 'LAST:%s'%chain2))
            aln.append_model(mdl_2, align_codes=alignCode, atom_files=alignCode+'.pdb')
            aln.align2d()
            aln.write(file='alignment.ali', alignment_format='PIR')

            """ select missing loop residues for modeling only"""
            class MyModel(automodel):
                # overide select_atom function, select only some residues
                def select_atoms(self):
                    # Select residues loopRange[0] to loopRange[1] (PDB numbering)
                    return selection(self.residue_range(str(loopRange[0])+":"+str(chain1),
                                                        str(loopRange[1])+":"+str(chain1)
                                                        )
                                     )

            a = MyModel(env,
                        alnfile='alignment.ali',
                        knowns=alignCode,
                        sequence=pdbCode,
                        )
            # a.auto_align()
            # get an automatic loop building
            a.make()

            # obtain successfully modeled pdb file names
            return (self.selectBestModeller(a.outputs))

    def addLoopsSimpleModel(self, pdbIn, chain1,
                            fastaSeq,
                            loopRanges,
                            no_lig=0,
                            verbose=False
                            ):
        '''
        add loop using no template
        :param pdbIn:
        :param chain1:
        :param fastaSeq:
        :param loopRanges:
        :param no_lig:
        :param verbose:
        :return:
        '''
        if pdbIn[-4:] == ".pdb" :
            pdbIn = pdbIn[:-4]

        env = environ()
        env.io.atom_files_directory = ['.', '../atom_files']

        if verbose :
            log.verbose()

        m = model(env, file=pdbIn)
        aln = alignment(env)
        aln.append_model(m, align_codes=pdbIn)
        aln.write(file=pdbIn + '.seq')

        ## add full sequence into the alignment file
        with open("alignment.ali",'wb') as tofile :
            with open(pdbIn + '.seq') as lines :
                tofile.write(lines)
            tofile.write('\nP1;%s_fill \n'%pdbIn)
            tofile.write('sequence::::::::: \n')
            for i in range(len(fastaSeq)) :
                tofile.write(fastaSeq[i])
                if i % 74 == 0 and i != 0 :
                    tofile.write("\n")

            tofile.write("%s"%"."*no_lig)
            tofile.write("*\n")

        class MyModel(automodel):
            def select_atoms(self):
                s = []
                for i in range(len(loopRanges)) :
                    s.append(selection(self.residue_range(
                        str(loopRanges[i][0]),
                        str(loopRanges[i][1])))
                    )
                return s

        a = MyModel(env, alnfile='alignment.ali', knowns = pdbIn, sequence = pdbIn+"_fill")
        a.starting_model = 1
        a.ending_model = 1
        #a.loop.md_level = refine.fast

        a.make()

        # obtain successfully modeled pdb file names
        return (self.selectBestModeller(a.outputs))

    def addLoopRefineModel(self, pdbIn, chain1,
                            fastaSeq,
                            loopRanges,
                            no_lig=0,
                            verbose=False):
        '''
        model the loops and refine them without any templates
        :param pdbIn:
        :param chain1:
        :param fastaSeq:
        :param loopRanges:
        :param no_lig:
        :param verbose:
        :return: the file name of the best model
        '''
        env = environ()
        if verbose :
            log.verbose()

        m = model(env, file=pdbIn)
        aln = alignment(env)
        aln.append_model(m, align_codes=pdbIn)
        aln.write(file=pdbIn + '.seq')

        ## add full sequence into the alignment file
        with open("alignment.ali", 'wb') as tofile:
            with open(pdbIn + '.seq') as lines:
                for s in lines :
                    if s[:-1] == "*" and s[0] != ">" :
                        tofile.write(s[:-1]+"."*no_lig+"*\n")
                    else :
                        tofile.write(s)
            tofile.write('\nP1;%s_fill \n' % pdbIn)
            tofile.write('sequence::::::::: \n')
            for i in range(len(fastaSeq)):
                tofile.write(fastaSeq[i])
                if i % 74 == 0 and i != 0:
                    tofile.write("\n")

            tofile.write("%s" % "." * no_lig)
            tofile.write("*\n")

        a = loopmodel(env, alnfile='alignment.ali', knowns = pdbIn , sequence = pdbIn+'_fill')

        a.starting_model = 1
        a.ending_model = 1

        a.loop.starting_model =  1
        a.loop.ending_model = 10
        a.loop.md_level = refine.fast

        a.make()

        # obtain successfully modeled pdb file names
        return self.selectBestModeller(a.outputs)

    def selectBestModeller(self, modellerOutput, verbose=False):
        # obtain successfully modeled pdb file names
        ok_models = [x for x in modellerOutput if x['failure'] is None]

        # Rank the models by DOPE score
        key = 'DOPE score'
        if sys.version_info[:2] == (2, 3):
            # Python 2.3's sort doesn't have a 'key' argument
            ok_models.sort(lambda a, b: cmp(a[key], b[key]))
        else:
            ok_models.sort(key=lambda a: a[key])

        # Get top model
        finalPDB = ok_models[0]
        if verbose:
            print("Top model: %s (DOPE score %.3f)" % (finalPDB['name'], finalPDB[key]))
        return (finalPDB)