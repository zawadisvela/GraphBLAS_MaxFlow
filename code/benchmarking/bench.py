
import sys
import os
import subprocess
from subprocess import PIPE

url_file_path = sys.argv[1]

file = open(url_file_path)

for url in file:


    output = subprocess.run(["ls"], stdout=PIPE, stderr=PIPE)
    print (output.stdout)
    break

    url = url.rstrip()
    command = 'wget '+url
    print(command)
    os.system(command)

    tar_name = url.split('/')[-1]
    command = 'tar -xf '+tar_name
    print(command)
    os.system(command)

    command = 'rm '+tar_name
    print(command)
    os.system(command)

    graph_name = tar_name.split('.')[0]
    filename = graph_name+'/'+graph_name+'.mtx'

    command = 'make edmund-karp'
    print(command)
    os.system(command)

    command = os.getcwd()+'/edmund-karp '+filename
    print(command)
    #os.system(command)

    command = 'rm -rf '+graph_name+'/'
    print(command)
    os.system(command)

    print()
    break
