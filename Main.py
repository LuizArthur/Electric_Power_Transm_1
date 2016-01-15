__author__ = 'LuizArthur'

import os
import time
import Calc
import matplotlib.pyplot as pp

tempo = time.clock()

main = os.getcwd()
Input = os.path.join(main, "Input")
Output = os.path.join(main, "Output")

if not os.path.exists(Output):
    os.mkdir(Output)

file_name = []
file_initiate_name = os.path.join(Input, "Input.txt")
file_initiate = open(file_initiate_name, "r")

for line in file_initiate:
    file_name.append(line.strip().split(" "))

file_initiate.close()

count_plot = 0

for name in file_name[0]:
    #print(name)
    path_name = os.path.join(Input, name)
    #print(path_name)
    if os.path.exists(path_name):
        with open(path_name) as file:
            Var_txt = []
            for line in file:
                Var_txt.append(line.strip().split(" "))

        file.close()

        calc = Calc.CALC(Var_txt, count_plot, name)
        count_plot = calc.show_count_plot()
    else:
        print("File "+path_name+" nao existe")

print("\nTempo de execucao: ", time.clock()-tempo)
pp.show()
