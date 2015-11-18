import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:",["ifile=","ofile=","n_phi="])
    except getopt.GetoptError:
        print 'test.py -i <input_run_name> -o <outputfile>'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            usage()
#             print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_run = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-p","--n_phi"):
            number_qs = arg
    print 'Input file is %s' % input_run
    print 'Output file is %s'% outputfile
    print 'Number of phis used is %s'%number_qs

if __name__ == "__main__":
   main(sys.argv[1:])
