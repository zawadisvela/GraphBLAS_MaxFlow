
import sys
import os
import subprocess
from subprocess import PIPE

thread_configs = [80, 40, 20] #, 8, 4, 1, 2]

runs = 3

def main():
    url_file_path = sys.argv[1]

    file = open(url_file_path)

    if not os.path.isfile('Makefile'):
        print('Run from build folder you dweeb')
        return

    if not os.path.isdir('output'):
        os.system('mkdir output')

    for url in file:

        #output = subprocess.run(["ls"], stdout=PIPE, stderr=PIPE)
        #print (output.stdout)
        #break

        command = 'cat '+url.rstrip()
        output = subprocess.run(command.split(), stdout=PIPE)
        print (output.stdout.decode("utf-8"))

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
        #filename = '../input/uk/uk.mtx'

        command = 'make edmund-karp'
        print(command)
        os.system(command)

        output_dir = 'output/'+graph_name
        if not os.path.isdir(output_dir):
            command = 'mkdir '+output_dir
            print(command)
            os.system(command)

        s_t_search_file = output_dir+'/'+graph_name+'s_t_set.txt'

        if True: #not os.path.isfile(s_t_search_file):
            command = 'make s_t_search'
            print(command)
            os.system(command)

            handle = open(s_t_search_file, 'w')
            command = './s_t_search '+filename
            print(command)
            output = subprocess.run(command.split(), stdout=handle)
            handle.close()

        s = -1
        t = -1
        s_t_file_handle = open(s_t_search_file, 'r')
        for line in s_t_file_handle:
            if line.startswith('s-t'):
                s,t = line.split(' ')[1].split('-')

        for threads in thread_configs:
            #command = 'export OMP_NUM_THREADS='+str(threads)
            #print(os.environ["OMP_NUM_THREADS"])
            os.environ["OMP_NUM_THREADS"] = str(threads)
            print(os.environ["OMP_NUM_THREADS"])
            #print(command)
            #os.system(command)
            #output = subprocess.run(command.split())

            out_name = output_dir+'/'+graph_name+'_results_'+str(threads)+'_threads.txt'
            if True: #not os.path.isfile(out_name):
                outfile = open(out_name, 'w')
                command = './edmund-karp '+filename+' '+str(runs)+' '+str(s)+' '+str(t)
                print(command)
                outfile.write(command+'\n'+'OMP_NUM_THREADS='+os.environ["OMP_NUM_THREADS"]+'\n\n')
                output = subprocess.run(command.split(), stdout=outfile)
                outfile.close()
                print('FINISHED RUN')
                #Adding average, max and min time
                outfile = open(out_name, 'r')
                times = []
                for line in outfile:
                    if line.startswith('Time'):
                        print('whaaat')
                        times.append(float(line.split(' ')[1]))
                    print(times)
                print('times:', times)
                outfile.close()
                print('times:', times)
                outfile = open(out_name, 'a')
                outfile.write('\nMin time: '+str(min(times)))
                outfile.write('\nMax time: '+str(max(times)))
                outfile.write('\nAvg time: '+str(sum(times)/len(times)))
                outfile.close()


        command = 'rm -rf '+graph_name+'/'
        print(command)
        os.system(command)

        print()


main()
