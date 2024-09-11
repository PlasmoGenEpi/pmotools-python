#!/usr/bin/env python3


import urllib.request, urllib.parse, urllib.error, os, shutil, tarfile, multiprocessing, subprocess, sys, socket
from pmotools.utils.color_text import ColorText as CT


class Utils:
    """
    A small utility class to hold static methods for various 1-off tools
    """

    @staticmethod
    def isMac():
        return sys.platform == "darwin"

    @staticmethod
    def connectedInternet():
        # from http://stackoverflow.com/questions/20913411/test-if-an-internet-connection-is-present-in-python
        try:
            # see if we can resolve the host name -- tells us if there is
            # a DNS listening
            host = socket.gethostbyname("www.google.com")
            # connect to the host -- tells us if the host is actually
            # reachable
            s = socket.create_connection((host, 80), 2)
            return True
        except:
            pass
        return False

    @staticmethod
    def which(program):
        # from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
        def is_exe(fnp):
            return os.path.isfile(fnp) and os.access(fnp, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
        return None

    @staticmethod
    def hasProgram(program):
        whichOutput = Utils.which(program)
        return None != whichOutput

    @staticmethod
    def run_in_dir(cmd, d):
        # print CT.boldBlack("here")
        cmd = "cd " + Utils.shellquote(d) + " && " + cmd + " && cd -"
        # print CT.boldBlack("newcmd")
        print(CT.boldGreen(cmd))
        Utils.run(cmd)

    @staticmethod
    def run(cmd):
        # from http://stackoverflow.com/a/4418193
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # output, errors = process.communicate()
        # sys.stdout.write(output.decode('utf-8'))
        # sys.stdout.flush()
        output = "";
        while True:
            nextline = process.stdout.readline().decode('utf-8')
            if nextline == '' and process.poll() != None:
                break
            sys.stdout.write(nextline)
            output = output + nextline
            sys.stdout.flush()
        exitCode = process.returncode
        if (exitCode == 0):
            return output
        raise Exception(cmd, exitCode, output)

    @staticmethod
    def runAndCapture(cmd):
        # from http://stackoverflow.com/a/4418193
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output, errors = process.communicate()
        # this is suppose to capture the output but it isn't for some reason so capturing it with the above
        exitCode = process.returncode
        if (exitCode == 0):
            return output.decode('utf-8')
        raise Exception(cmd, exitCode, output.decode('utf-8'), errors)

    @staticmethod
    def shellquote(s):
        # from http://stackoverflow.com/a/35857
        return "'" + s.replace("'", "'\\''") + "'"

    @staticmethod
    def num_cores():
        return multiprocessing.cpu_count()

    @staticmethod
    def mkdir(d):
        """mkdir if it doesn't already exist """
        if not os.path.exists(d):
            print(CT.boldText("mkdir"), CT.boldGreen(d))
            os.makedirs(d)

    @staticmethod
    def get_file(url, d):
        """get file from url and put it into directory d, return new name """
        fn = url.split('/')[-1]
        out_fnp = os.path.join(d, fn)
        urllib.request.urlretrieve(url, out_fnp)
        return out_fnp

    @staticmethod
    def get_file_if_size_diff(url, d):
        """only download the file if it's needed, not completely fail proof since it is
        just a size check but fairly likely not to be the same for a difference """
        fn = url.split('/')[-1]
        out_fnp = os.path.join(d, fn)
        net_file_size = int(urllib.request.urlopen(url).info()['Content-Length'])
        if os.path.exists(out_fnp):
            fn_size = os.path.getsize(out_fnp)
            if fn_size == net_file_size:
                print("skipping download of", CT.boldGreen(fn))
                return out_fnp
            else:
                print("files sizes differed:", "on disk:", fn_size, "from net:", net_file_size)
        print("retrieving", CT.boldGreen(fn), "from", CT.boldBlue(url))
        urllib.request.urlretrieve(url, out_fnp)
        return out_fnp

    @staticmethod
    def rm_rf(d):
        """remove directory forcibly"""
        if os.path.exists(d):
            print(CT.boldText("rm -rf"), CT.boldRed(d))
            shutil.rmtree(d)

    @staticmethod
    def untar(fnp, d):
        """ un pack compressed file, guessing format based on extension"""
        if fnp.endswith(".tar.gz"):
            tar = tarfile.open(fnp, "r:gz")
        elif fnp.endswith(".tgz"):
            tar = tarfile.open(fnp, "r:gz")
        elif fnp.endswith(".tar.bz2"):
            tar = tarfile.open(fnp, "r:bz2")
        elif fnp.endswith(".tar"):
            tar = tarfile.open(fnp, "r")
        else:
            raise Exception("invalid file? " + fnp)
        print("untarring", CT.boldGreen(fnp), "to", CT.boldBlue(d))
        tar.extractall(d)
        tar.close()

    @staticmethod
    def getStrFromStrOrList(inputArg):
        if isinstance(inputArg, list):
            return str(inputArg[0])
        elif not isinstance(inputArg, str):
            return str(inputArg)
        else:
            return inputArg

    @staticmethod
    def clear_dir(d):
        """ forcibly delete directory and then re-make it"""
        Utils.rm_rf(d)
        Utils.mkdir(d)

    @staticmethod
    def appendStrAsNeeded(input : str, ending : str):
        """
        if a string doesn't end with a specific ending, append it, this is useful for ensuring file extensions are in output names without accidentally doubling it

        :param input: the string to be appended
        :param ending: the desired ending
        :return: the string with eh ending appended if it doesn't already end with it
        """
        if not input.endswith(ending):
            return input + ending
        return input

    @staticmethod
    def parse_delimited_input_or_file(input : str, delim : str = ",") -> list[str]:
        """
        If the input is a file name then read in each line of the file for the argument, otherwise return a list of items delimited by delimiter

        :param input: the argument to parse or a name of a file to read in
        :param delim: the delimiter to split on
        :return: a list of strings
        """
        ret = []
        if len(input) <=255 and os.path.exists(input):
            with open(input) as file:
                ret = [line.rstrip() for line in file]
        else:
            ret = input.split(delim)
        return ret

    @staticmethod
    def process_delimiter_and_output_extension(delim : str, output_extension : str = ".txt"):
        """
        Process delimiter and extension, this allows for delim to be listed as tab or comma and it will replace appropriately the

        :param delim: the delimiter to process
        :param output_extension: the output extension
        :return: delimiter, extension
        """
        out_delim = delim
        out_output_extension = output_extension

        output_extension = ".txt"
        if delim == "tab" or delim == "\t":
            out_delim = "\t"
            out_output_extension = ".tsv"
        elif delim == "comma" or delim == ",":
            out_delim = ","
            out_output_extension = ".csv"
        return out_delim, out_output_extension

    @staticmethod
    def outputfile_check(output_file : str, overwrite : bool = False):
        """
        Check to see if the output file exists if overwrite is turned on or not

        :param output_file: the output file that will be written to
        :param overwrite: whether or not the output file can be overwritten
        :return: None
        """
        # only overwrite an existing file if --overwrite is on
        if "STDOUT" != output_file and os.path.exists(output_file) and not overwrite:
            raise Exception(
                "Output file " + output_file + " already exists, use --overwrite to overwrite it")

    @staticmethod
    def inputOutputFileCheck(inputFile : str, outputFile : str, overwrite : bool = False):
        """
        Check to see if an input file exists and if the output file exists if overwrite is turned on or not

        :param inputFile: the file that will be read in
        :param outputFile: the output file that will be written to
        :param overwrite: whether or not the output file can be overwritten
        :return: None
        """
        # make sure file exists
        if not os.path.exists(inputFile):
            raise FileNotFoundError(inputFile)

        Utils.outputfile_check(outputFile, overwrite)

    @staticmethod
    def inputOutputFileCheckFromArgParse(args):
        """
        Check if input file exists and if output file exist check if --overwrite flag is set

        :param args: parsed arguments from argparse, expecting the fields 'file', 'output', and 'overwrite'
        :return:
        """
        Utils.inputOutputFileCheck(args.file, args.output, args.overwrite)
