
import sys
import os
import subprocess
from subprocess import PIPE

thread_configs = [80, 40, 20] #, 8, 4, 1, 2]

runs = 5

def main():
    url_file_path = sys.argv[1]

    file = open(url_file_path)

    if not os.path.isfile('Makefile'):
        print('Run from build folder you dweeb')
        return

    if not os.path.isdir('output'):
        os.system('mkdir output')

    for url in file:

        #output = subprocess.run(['ls'], stdout=PIPE, stderr=PIPE)
        #print (output.stdout)
        #break

        command = 'cat '+url.rstrip()
        output = subprocess.run(command.split(), stdout=PIPE)
        print (output.stdout.decode('utf-8'))

        input_dir = '../'
        graph_name = url.split('/')[-1].split('.')[0]
        filename = input_dir+graph_name+'/'+graph_name+'.mtx'

        output_dir = 'output/'+graph_name
        if not os.path.isdir(output_dir):
            command = 'mkdir '+output_dir
            print(command)
            os.system(command)

        all_done = True
        for threads in thread_configs:
            out_name = output_dir+'/'+graph_name+'_results_'+str(threads)+'_threads.txt'
            if not os.path.isfile(out_name) or len(open(out_name).readlines()) == 0:
                all_done = False

        if all_done:
            continue

        print("Looking for", filename)

        if not os.path.isfile(filename):
            print("Didn't find", filename)
            url = url.rstrip()
            command = 'wget '+url
            print(command)
            os.system(command)

            tar_name = url.split('/')[-1]
            command = 'tar -C '+input_dir+' -xf '+tar_name
            print(command)
            os.system(command)

            command = 'rm '+tar_name
            print(command)
            os.system(command)

        #filename = '../input/uk/uk.mtx'

        command = 'make edmund-karp'
        print(command)
        os.system(command)

        s_t_search_file = output_dir+'/'+graph_name+'s_t_set.txt'

        if not os.path.isfile(s_t_search_file) or len(open(s_t_search_file).readlines()) == 0:
            command = 'make s_t_search'
            print(command)
            os.system(command)

            handle = open(s_t_search_file, 'w')
            command = './s_t_search '+filename
            print(command)
            os.environ['OMP_NUM_THREADS'] = str(80)
            #output = subprocess.run(command.split(), stdout=handle)
            os.system(command)
            handle.close()

        s = -1
        t = -1
        s_t_file_handle = open(s_t_search_file, 'r')
        for line in s_t_file_handle:
            if line.strip().startswith('s-t'):
                s,t = line.split(' ')[1].split('-')

        for threads in thread_configs:
            #command = 'export OMP_NUM_THREADS='+str(threads)
            #print(os.environ['OMP_NUM_THREADS'])
            os.environ['OMP_NUM_THREADS'] = str(threads)
            print(os.environ['OMP_NUM_THREADS'])
            #print(command)
            #os.system(command)
            #output = subprocess.run(command.split())

            out_name = output_dir+'/'+graph_name+'_results_'+str(threads)+'_threads.txt'
            if not os.path.isfile(out_name) or len(open(out_name).readlines()) < 10:
                if os.path.isfile(out_name):
                    print('File was only:', len(open(out_name).readlines()), 'long')
                    command = 'mv '+out_name+' '+out_name.split('.')[0]+'_ERROR.txt'
                    print(command)
                    subprocess.run(command.split())

                outfile = open(out_name, 'w')
                command = './edmund-karp '+filename+' '+str(runs)+' '+str(s)+' '+str(t)
                print(command)
                outfile.write('\n'+command+'\n'+'OMP_NUM_THREADS='+os.environ['OMP_NUM_THREADS']+'\n\n')
                output = subprocess.run(command.split(), stdout=outfile, stderr=outfile)
                """
                p = subprocess.Popen(command.split(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT)
                for line in iter(p.stdout.readline, b''):
                    outfile.write(line.decode("utf-8"))
                    print(">>> " + line.decode("utf-8").rstrip())
                """
                outfile.close()
                print('FINISHED RUN')
                #Adding average, max and min time
                outfile = open(out_name, 'r')
                times = []
                for line in outfile:
                    if line.startswith('Total run time'):
                        times.append(float(line.split(':')[1]))
                #print('times:', times)
                outfile.close()
                print('times:', times)
                outfile = open(out_name, 'a')
                try:
                    outfile.write('\nTimes: '+str(times))
                    outfile.write('\nMin time: '+str(min(times)))
                    outfile.write('\nMax time: '+str(max(times)))
                    outfile.write('\nAvg time: '+str(sum(times)/len(times)))
                except:
                    outfile.write('Error while reading times')
                outfile.close()


        #command = 'rm -rf '+input_dir+graph_name+'/'
        #print(command)
        #os.system(command)
        print("NOT deleting the graph, got disk to spare")

        print()


main()
