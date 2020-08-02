import linecache
from collections import OrderedDict

class TopModifier:
    def __init__(self, top_in):
        self.topin = top_in

    def addFromLines(self, top_out, lines, field, after=True):
        '''
        input the topol file, and add some lines after a field
        :param top_out: str, output file name
        :param lines: list of str, the lines to be added
        :param after_field: str, field name where new lines to be added
        :return: 1
        '''

        nlines = self.getFileLineNums(self.topin)

        with open(top_out, 'w') as tofile:
            field_linenum = self.topFields()[field]

            if after:
                head = field_linenum + 1
                tail = field_linenum + 1
            else:
                head = field_linenum - 1
                tail = field_linenum

            for i in range(head):
                tofile.write(linecache.getline(self.topin, i))

            # for line in lines :
            #    tofile.write(line)
            tmp = [tofile.write(x) for x in lines]

            for j in range(tail, nlines + 1):
                tofile.write(linecache.getline(self.topin, j))

        return 1

    def addFromFile(self, top_out, newline_file, field, after=True):
        '''
        add the lines in a file to the topol file
        :param top_out: str, output file name
        :param newline_file: str, the file where the lines would be added to newfile
        :param after_field: str, the field after which new lines will be added
        :return:
        '''

        try:
            with open(newline_file) as lines:
                newlines = lines.readlines()
        except FileNotFoundError:
            newlines = ""

        return self.addFromLines(top_out, newlines, field, after)

    def getFileLineNums(self, filein):
        '''
        get number of lines from a given file
        :param filein: str, input file name
        :return:
        '''

        try:
            with open(filein) as lines:
                nlines = len(lines.readlines())
            return nlines

        except FileNotFoundError:
            return 0

    def topFields(self):
        '''
        input a topol file, or a gromacs itp file,
        return the fields and their line numbers
        use sorted dictionary
        :return:
        '''

        fields, line_no = [], []
        nlines = int(self.getFileLineNums(self.topin))

        for i in range(nlines):

            line = linecache.getline(self.topin, i)

            if "[" in line and "]" in line and line.split()[0][0] != ";":
                fields.append(line.split()[1])
                line_no.append(i)
            else:
                pass

        fields_lineno = OrderedDict()

        for k in range(len(fields)):
            fields_lineno[fields[k]] = line_no[k]

        # print(fields)
        return fields_lineno


class IndexModifier(TopModifier):
    def appendFields(self, fields=[], output="output.ndx", field_name="Combined"):
        '''
        cat several fields together at the end of a output file
        :param fields: list of str, the field names
        :param output:
        :param field_name:
        :return:
        '''

        all_lines = []

        for f in fields:
            flines = self.getFieldContents(f)

            all_lines += flines

        if len(field_name):
            fn = field_name
        else:
            fn = "_".join(fields)

        all_lines = ["[ {:s} ] \n".format(fn), ] + all_lines

        all_lines = [x for x in all_lines if len(x.split())]

        return self.addLinesFileEnd(output, all_lines)

    def getFieldContents(self, field):
        '''
        get the contents (lines) of a field
        :param field: str, field name
        :return: list of str, lines
        '''
        fields_nlines = self.topFields()
        all_fields = list(fields_nlines.keys())

        thisfield_ln = fields_nlines[field]
        nextfield_ln = fields_nlines[all_fields[all_fields.index(field) + 1]]

        field_content = []
        for ln in range(thisfield_ln + 1, nextfield_ln):
            field_content += [linecache.getline(self.topin, ln)]

        return field_content

    def addLinesFileEnd(self, output, lines):
        '''
        append lines at the end of a file
        :param output: str, output file
        :param lines: list of str
        :return:
        '''

        with open(output, 'w') as tofile:
            prev_lines = open(self.topin, 'r').readlines()
            tmp = [tofile.write(x) for x in prev_lines]

            tmp = [tofile.write(x) for x in lines]

        tofile.close()

        return 1
