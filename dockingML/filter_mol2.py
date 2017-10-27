#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .mol2IO import Mol2IO
import os

class FilterMol2(Mol2IO) :

    def extractInfor(self, input, property, lineindex=1):
        """
        from the mol2 file extract necessary information
        :param input:
        :param property:
        :param lineindex:
        :return:
        """
        with open(input) as lines :
            condition = False
            count     = -1
            for s in lines :
                if s.strip("\n") == property :
                    condition = True

                if condition :
                    count += 1
                    if count == lineindex :
                        return s.strip("\n")
        return "9999"

    def compareProperty(self, input, property, cutoff, largerthan=True, valueIndex=0):
        """
        compare the property in a mol2 file to the cutoff
        :param input:
        :param property:
        :param cutoff:
        :param largerthan:
        :param valueIndex:
        :return:
        """
        #value = self.triposInfor(input=input, property=property)[valueIndex]
        value  = self.extractInfor(input, property, valueIndex)

        def toFloat(x) :
            newx = 9999.0
            if ">" in x or "<" in x :
                newx = str(x.strip(">"))
                newx = str(newx.strip("<"))
            else :
                newx = str(x)

            if "," in newx :
                raw_v = newx.split(",")
                newx = 0
                for i in raw_v :
                    newx += (1000 ** (len(raw_v) - i - 1)) * float(raw_v[i])

            else :
                if newx :
                    newx = float(newx)
                else :
                    newx = 9999.9

        value = toFloat(value)

        if largerthan :
            if value >=  cutoff:
                return True
            else : return False
        else :
            if value <= cutoff :
                return True
            else : return False

    def catListOfFiles(self, input_list, output):
        """
        cat a list of files to a new file
        :param input_list:
        :param output:
        :return:
        """

        os.system("touch temp")

        for filen in input_list :
            self.catenateMol("temp", filen, "temp")

        os.rename("temp", output)
        #os.remove("temp")

