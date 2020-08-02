
from mdanaly import extract
import subprocess as sp
from automd.utils import gentop


class CleanPDB :
    '''
    Prepare the pdb file, add missing atoms, add hydrogens, add missing loops
    '''
    def __init__(self, filename, obabel='obabel'):
        self.pdb = filename
        # the input file is not a pdb, or pdbqt file, use obabel convert it
        if filename.split(".")[-1] not in ["pdb","pdbqt"] :
            gpdb = gentop.GenerateTop()
            gpdb.runObabel(obabelexe=obabel, input=filename, output=filename+".pdb")
            self.pdb = filename+".pdb"

    def extractFrame(self, frame=1, headersave=True) :
        '''
        extract a frame (default is first) of the pdb file

        :param frame:
        :param headersave:
        :return: the specific frame from a large pdb file
        '''
        extractPDB = extract.ExtractPDB()
        extractPDB.extract_frame(self.pdb, self.pdb[:-5]+"_%d.pdb"%frame, no_frames=[frame])
        return(self.pdb[:-5]+"_%d.pdb"%frame)

    def processHETATM(self, filename,
                      hetatmsave=['WAT','HOH'],
                      dropWater=True,
                      cleanedPDB='cleaned_', headersave=True,
                      selectedChains = [],
                     ):
        '''
        :param filename:
        :param hetatmsave:
        :param dropWater: remove water molecules
        :param cleanedPDB:
        :param headersave: save the pdb header information
        :param selectedChains:
        :return: the pdbfile after removing unnecessary information
        '''
        tofile = open(cleanedPDB+filename, 'wb')
        with open(filename) as lines :
            for s in lines :
                if len(s.split()) > 0 and \
                                s.split()[0] in ['ATOM', 'HETATM', 'TER', 'END', 'ENDMDL'] :
                    if dropWater :
                        if s[21] in selectedChains and s[17:20].strip() not in ['HOH', 'WAT'] \
                                and s[17:20].strip() in hetatmsave :
                            tofile.write(s)
                    else :
                        if s[21] in selectedChains and s[17:20].strip() in hetatmsave :
                            tofile.write(s)
                else :
                    if headersave :
                        tofile.write(s)

        return(cleanedPDB+filename)

    def removeLonePair(self, inpdb, outpdb):
        '''
        input a pdb file, remove the lone pair electron lines in the file
        :param inpdb:
        :param outpdb:
        :return:
        '''

        job = sp.Popen("awk \'$1 ~ /AT/ && $NF !~ /Xx/ {print $0}\' %s > %s "%(inpdb, outpdb), shell=True)
        job.communicate()

        return 1
