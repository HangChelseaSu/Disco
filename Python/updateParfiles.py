import sys
import argparse as ap

def loadTemplate(filename):
    f = open(filename, "r")

    template = []
    for line in f:
        if line == "" or line.isspace():
            continue
        words = line.split()
        if len(words) >= 2 and words[0][0].isalpha():
            # Assuming this is a parameter
            par = words[0]
            template.append((par,line))
        else:
            template.append((line,))

    f.close()

    return template

def updateParfile(parfile, template):

    f = open(parfile, "r")
    lines = f.readlines()
    f.close()

    n = 0

    count0 = len(lines)

    for tline in template:
        if len(tline) == 1:
            if tline[0] not in lines:
                lines.insert(n, tline[0])
                n += 1
            else:
                n = lines.index(tline[0])+1
        else:
            found = False
            for i, line in enumerate(lines):
                words = line.split()
                if len(words) < 2:
                    continue
                if words[0] == tline[0]:
                    found = True
                    n = i+1
                    break
            if found is False:
                lines.insert(n, tline[1])
                n += 1

    for i in range(len(lines)):
        if lines[i][-1] != '\n':
            lines[i] += '\n'

    count1 = len(lines)
    print("Added {0:d} lines".format(count1-count0))

    f = open(parfile, "w")
    f.writelines(lines)
    f.close()



if __name__ == "__main__":

   parser = ap.ArgumentParser(description="Updates parameter files to include parameters from template.")
   parser.add_argument('template', help="Template parameter file.")
   parser.add_argument('parfiles', nargs='+', help="Parameter files to update. Every parameter in the template not present in these files will be written in to the files.")

   args = parser.parse_args()

   template = args.template
   parfiles = args.parfiles

   template_str = loadTemplate(template)

   for file in parfiles:
       updateParfile(file, template_str)
