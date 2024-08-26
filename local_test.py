import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--num_parties',  action='store', dest='n', type = int, required = True)
parser.add_argument('-c', '--command',  action='store', dest='c', type = str, required = True)
parser.add_argument('-s', '--start_party',  action='store', dest='s', type = int, required = True)
parser.add_argument('-e', '--end_party',  action='store', dest='e', type = int, required = True)
parser.add_argument('-p', '--num_processes',  action='store', dest='p', type = int, required = False)
args = parser.parse_args()

def test(n):
    s = ""
    for i in range(args.s, args.e):
        # if i == 5: continue
        if args.c == 'yao':
            s +=  "./bin/test_"+ args.c +" " + ( i + 1).__str__() + " 12345 test/and.txt  " + n.__str__() + " ip_file.txt"
        if args.c == 'mp_circuit':
            s +=  "./bin/test_"+ args.c +" " + ( i + 1).__str__() + " 12345 ../emp-aby/modsum.txt  " + n.__str__()
        else:
            s += "./bin/test_"+ args.c +" " + ( i + 1).__str__() + " 12345 " + n.__str__() + " ip_file.txt"
        if i != args.e-1: 
            s += " & "
    print(s)
    return s

sleep_time_batch = 70
#test(args.n)
if args.p:
    s = ""
    sleep_time = 0
    for i in range(args.p):
        s += "sleep " + (sleep_time * sleep_time_batch).__str__() + "; " + test(args.n)
        sleep_time += 1
        if i != args.p-1:
            s += " & "
        
    print(s)
    subprocess.Popen(s, shell="True")
else:
    subprocess.Popen(test(args.n), shell="True")
