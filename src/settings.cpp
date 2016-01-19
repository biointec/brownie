/***************************************************************************
 *   Copyright (C) 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)     *
 *   Copyright (C) 2014, 2015 Mahdi Heydari (mahdi.heydari@intec.ugent.be) *
 *   This file is part of Brownie                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "settings.h"
#include "library.h"

#include <iostream>
#include <fstream>
//#include <libxml/parser.h>
using namespace std;

// ============================================================================
// SETTINGS CLASS PRIVATE
// ============================================================================

void Settings::printProgramInfo() const
{
        cout << "Brownie read error correction v." << BROWNIE_MAJOR_VERSION << "."
             << BROWNIE_MINOR_VERSION << "." << BROWNIE_PATCH_LEVEL << "\n";

        cout << "Copyright 2014, 2015 Jan Fostier (jan.fostier@intec.ugent.be)\n";
        cout << "Copyright 2014, 2015 Mahdi Heydari (mahdi.heydari@intec.ugent.be)\n";
        cout << "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE.\n" << endl;

        cout << "Compilation settings:\n";
        cout << "  MAXKMERLENGTH = " << MAXKMERLENGTH << "\n";
#ifdef HAVE_ZLIB
        cout << "  ZLIB support = enabled";
#else
        cout << "  ZLIB support = disabled";
#endif
        cout << endl;
}

void Settings::printUsage() const
{
        cout << "Usage: brownie command [options] [file_options] file1 [[file_options] file2]...\n";
        cout << "Corrects sequence reads in file(s)\n\n";

        cout << " [options]\n";
        cout << "  -h\t--help\t\t\tdisplay help page\n";
        cout << "  -i\t--info\t\t\tdisplay information page\n";
        cout << "  -s\t--singlestranded\tenable single stranded DNA [default = false]\n\n";

        cout << "  command\n";
        cout << "  Assembly\t\t\tDe novo assembly\n";
        cout << "  ErrorCorrection\t\tDe novo error correction\n";
        cout << "  GraphConstruction\t\tBuild de bruijn graph\n";
        cout << "  GraphCorrection\t\tCorrect constructed de bruijn graph\n\n";


        cout << " [options arg]\n";
        cout << "  -k\t--kmersize\t\tkmer size [default = 31]\n";
        cout << "  -t\t--threads\t\tnumber of threads [default = available cores]\n";
        cout << "  -g\t--genomesize\t\tsize of the genome [default = auto]\n";
        cout << "  -v\t--visits\t\tmaximal number of visited nodes in one bubble detection [default = 1000]\n";
        cout << "  -d\t--depth\t\t\tmaximal number of visited nodes in one read correction [default = 1000]\n";
        cout << "  -e\t--essa\t\t\tsparseness factor of the enhanced sparse suffix array [default = 1]\n";
        cout << "  -c\t--cutoff\t\tvalue to separate ture and false nodes based on their coverage [default = calculated based on poisson mixture model]\n";

        cout << "  -p\t--pathtotmp\t\tpath to directory to store temporary files [default = current directory]\n\n";

        cout << " [file_options]\n";
        cout << "  -o\t--output\t\toutput file name [default = inputfile.corr]\n";

        cout << " examples:\n";
        cout << "  ./brownie inputA.fastq\n";
        cout << "  ./brownie -k 29 -t 4 -g 2800000 -o outputA.fasta inputA.fasta -o outputB.fasta inputB.fastq\n";
}

// ============================================================================
// SETTINGS CLASS PUBLIC
// ============================================================================

Settings::Settings() : kmerSize(31), numThreads(std::thread::hardware_concurrency()),
        genomeSize(0), doubleStranded(true), essa_factor(1), max_visits(1000),
        max_depth(1000),cutoff(0), skip_stage_4(false), skip_stage_5(false) {}

void Settings::parseCommandLineArgumentsEC(int argc, char** args,
                                         LibraryContainer& libCont)
{

        // parse all input arguments
        string inputFilename, outputFilename;
        for (int i = 2; i < argc; i++) {
                string arg(args[i]);

                if (arg.empty())
                        continue;              // this shouldn't happen

                 if ((arg == "-k") || (arg == "--kmersize")) {
                        i++;
                        if (i < argc)
                                kmerSize = atoi(args[i]);
                } else if ((arg == "-t") || (arg == "--threads")) {
                        i++;
                        if (i < argc)
                                numThreads = atoi(args[i]);
                } else if ((arg == "-g") || (arg == "--genomesize")) {
                        i++;
                        if (i < argc)
                                genomeSize = atoi(args[i]);
                } else if ((arg == "-e") || (arg == "--essa")) {
                        i++;
                        if (i < argc)
                                essa_factor = atoi(args[i]);
                } else if ((arg == "-v") || (arg == "--visits")) {
                        i++;
                        if (i < argc)
                                max_visits = atoi(args[i]);
                } else if ((arg == "-d") || (arg == "--depth")) {
                        i++;
                        if (i < argc)
                                max_depth = atoi(args[i]);
                } else if ((arg == "-c") || (arg == "--cutoff")) {
                        i++;
                        if (i < argc)
                                cutoff = atoi(args[i]);
                } else if ((arg == "-s") || (arg == "--singlestranded")) {
                        doubleStranded = false;
                } else if ((arg == "-p") || (arg == "--pathtotmp")) {
                        i++;
                        if (i < argc)
                                pathtotemp = args[i];
                } else if ((arg == "-o") || (arg == "--output")) {
                        i++;
                        if (i < argc)
                                outputFilename = args[i];
                } else {        // it must be an input file
                        inputFilename = args[i];
                        ReadLibrary lib = ReadLibrary(inputFilename, outputFilename);
                        libCont.insert(lib);
                        outputFilename.clear();
                }
        }
         if (!pathtotemp.empty()) {
                if ((pathtotemp.back() != '/') && (pathtotemp.back() != '\\'))
                        pathtotemp.push_back('/');
        }
           // final check: see if we can write to the temporary directory
        ofstream ofs(pathtotemp + "log.txt");
        if (!ofs.good()) {
                cerr << "brownie: cannot write to directory: " << pathtotemp << "\n";
                cerr << "Please make sure the path exists" << endl;
                exit(EXIT_FAILURE);
        }
        for (int i = 0; i < argc-1; i++)
                ofs << argc << " ";
        ofs << args[argc-1] << endl;
        ofs.close();
}
void Settings::parseCommandLineArgumentsMain(int argc, char** args,
                                             LibraryContainer& libCont)
{
        string arg(args[1]);

        if (arg== "ErrorCorrection"){
                parseCommandLineArgumentsEC(argc,args,libCont);
        }else if (arg== "Assembly"){

        }
        else if (arg== "GraphCorrection") {
                skip_stage_5 = true;
                parseCommandLineArgumentsEC(argc,args,libCont);
        }
        else if(arg=="GraphConstruction"){
                skip_stage_4 = true;
                skip_stage_5 = true;
                parseCommandLineArgumentsEC(argc,args,libCont);
        }
        else if ((arg == "-h") || (arg == "--help")) {
                printUsage();
                exit(EXIT_SUCCESS);
        } else if ((arg == "-i") || (arg == "--info")) {
                printProgramInfo();
                exit(EXIT_SUCCESS);
        }
        checkInputArguments(libCont);
}
void Settings::checkInputArguments(LibraryContainer& libCont){
                // perform sanity check on input parameters
        if (numThreads == 0)
                numThreads = 1;

        if (numThreads > std::thread::hardware_concurrency()) {
                cerr << "WARNING: number of threads is bigger than available number of cores" << endl;
        }

        if (kmerSize <= KMERBYTEREDUCTION *4) {
                cerr << "The kmer size must be at least " << 4*KMERBYTEREDUCTION + 1 << endl;
                throw ("Invalid argument");
        }

        if (kmerSize % 2 == 0) {
                cerr << "The kmer size must be odd" << endl;
                throw ("Invalid argument");
        }
        if (kmerSize > MAXKMERLENGTH) {
                size_t maxKmerLength = MAXKMERLENGTH;
                if (maxKmerLength % 2 == 0)
                        maxKmerLength--;
                cerr << "The kmer size can be at most " << maxKmerLength << endl;
                cerr << "Recompile Brownie with a higher MAXKMERLENGTH if a higher kmer size is desired" << endl;
                throw ("Invalid argument");
        }
        if (libCont.getSize() == 0) {
                cerr << "brownie: missing input read file\n";
                cerr << "Try 'brownie --help' for more information" << endl;
                exit(EXIT_FAILURE);
        }
}
