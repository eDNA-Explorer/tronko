#include "options.h"
#include <getopt.h>

static struct option long_options[]=
{
	{"help", no_argument, 0, 'h'},
	{"paired", no_argument, 0, 'p'},
	{"single", no_argument, 0, 's'},
	{"reference", no_argument, 0, 'r'},
	{"partition-directory", no_argument, 0, 'y'},
	{"reference-file", required_argument, 0, 'f'},
	{"tree-file", required_argument, 0, 't'},
	{"msa-file", required_argument, 0, 'm'},
	{"tax-file", required_argument, 0, 'x'},
	{"partitions-directory", required_argument, 0, 'd'},
	{"results", required_argument, 0, 'o'},
	{"single-read-file", required_argument, 0, 'g'},
	{"paired-read-file1", required_argument, 0, '1'},
	{"paired-read-file2", required_argument, 0, '2'},
	{"fasta-file", required_argument, 0, 'a'},
	{"Cinterval", required_argument, 0, 'c'},
	{"reverse-single-read",no_argument, 0, 'v'},
	{"reverse-paired-read",no_argument, 0, 'z'},
	{"number-of-cores",required_argument, 0,'C'},
	{"number-of-lines-to-read",required_argument, 0, 'L'},
	{"print-alignments",required_argument, 0, 'P'},
	{"use-nw",no_argument, 0, 'w'},
	{"print-unassigned",no_argument, 0, 'U'},
	{"use-leaf-portion",no_argument, 0, 'e'},
	{"padding",required_argument, 0, 'n'},
	{"fastq",no_argument, 0, 'q'},
	{"print-node-info",required_argument, 0, '5'},
	{"skip-bwa-build",no_argument,0, '6'},
	{"score-constant",required_argument,0, 'u'},
	{"print-all-scores",no_argument,0,'7'},
	{"print-alignments-dir",required_argument,0,'3'},
	{"tree-dir",required_argument,0,'4'},
	{"verbose",optional_argument,0,'V'},
	{"log-file",required_argument,0,'l'},
	{"enable-resource-monitoring",no_argument,0,'R'},
	{"enable-timing",no_argument,0,'T'},
	{"tsv-log",required_argument,0,0},  // Long option only, no short form
	{"early-termination",no_argument,0,0},
	{"no-early-termination",no_argument,0,0},
	{"strike-box",required_argument,0,0},
	{"max-strikes",required_argument,0,0},
	{"enable-pruning",no_argument,0,0},
	{"disable-pruning",no_argument,0,0},
	{"pruning-factor",required_argument,0,0},
	{0, 0, 0, 0}  // Terminating entry required by getopt_long
};

char usage[] = "\ntronko-assign [OPTIONS] -r -f [TRONKO-BUILD DB FILE] -a [REF FASTA FILE] -o [OUTPUT FILE]\n\
	\n\
	-h, usage:\n\
	-r, REQUIRED, use a reference\n\
	-f [FILE], REQUIRED, path to reference database file, can be gzipped\n\
	-a [FILE], REQUIRED, path to reference fasta file (for bwa database)\n\
	-o [FILE], REQUIRED, path to output file\n\
	-p, use paired-end reads\n\
	-s, use single reads\n\
	-v, when using single reads, reverse-complement it\n\
	-z, when using paired-end reads,  reverse-complement the second read\n\
	-g [FILE], compatible only with -s, path to single-end reads file\n\
	-1 [FILE], compatible only with -p, path to paired-end forward read file\n\
	-2 [FILE], compatible only with -p, path to paired-end reverse read file\n\
	-c [INT], LCA cut-off to use [default:5]\n\
	-C [INT], number of cores [default:1]\n\
	-L [INT], number of lines to read for assignment [default:50000]\n\
	-P, print alignments to stdout\n\
	-w, use Needleman-Wunsch Alignment Algorithm (default: WFA)\n\
	-q, Query is FASTQ [default is FASTA]\n\
	-e, Use only a portion of the reference sequences\n\
	-n [INT], compatible only with -e, Padding (Number of bases) to use in the portion of the reference sequences\n\
	-5 [FILE], Print tree number and leaf number and exit\n\
	-6, Skip the bwa build if database already exists\n\
	-u, Score constant [default: 0.01]\n\
	-7, Print scores for all nodes [scores_all_nodes.txt]\n\
	-V [LEVEL], Enable verbose logging [0=ERROR, 1=WARN, 2=INFO, 3=DEBUG] [default: disabled]\n\
	-l [FILE], Log file path [default: stderr only]\n\
	-R, Enable resource monitoring (memory/CPU usage)\n\
	-T, Enable timing information\n\
	--tsv-log [FILE], Export memory stats to TSV file for analysis\n\
	\n\
	Optimization Options:\n\
	--early-termination, Enable early termination during tree traversal\n\
	--no-early-termination, Disable early termination (default)\n\
	--strike-box [FLOAT], Strike box size as multiplier of Cinterval [default: 1.0]\n\
	--max-strikes [INT], Maximum strikes before termination [default: 6]\n\
	--enable-pruning, Enable subtree pruning\n\
	--disable-pruning, Disable subtree pruning (default)\n\
	--pruning-factor [FLOAT], Pruning threshold = factor * Cinterval [default: 2.0]\n\
	\n";

void print_help_statement(){
	printf("%s", &usage[0]);
	return;
}

void parse_options(int argc, char **argv, Options *opt){
	int option_index, success;
	int c;  // getopt_long returns int, not char
	
	if (argc==1){
		print_help_statement();
		exit(0);
	}
	while(1){
		c=getopt_long(argc,argv,"hpsrqw6yevUzP:75:f:u:t:m:d:o:x:g:1:2:a:c:n:3:4:C:L:V::l:RT",long_options, &option_index);
		
		// Handle end of options and errors
		if (c == -1) {
			break;
		}
		if (c == '?' || c == 255) {
			fprintf(stderr, "Unknown option or missing argument\n");
			print_help_statement();
			exit(1);
		}
		
		switch(c){
			case 0:
				// Handle long options without short equivalents
				if (strcmp(long_options[option_index].name, "tsv-log") == 0) {
					strncpy(opt->tsv_log_file, optarg, sizeof(opt->tsv_log_file) - 1);
					opt->tsv_log_file[sizeof(opt->tsv_log_file) - 1] = '\0';
				}
				else if (strcmp(long_options[option_index].name, "early-termination") == 0) {
					opt->early_termination = 1;
				}
				else if (strcmp(long_options[option_index].name, "no-early-termination") == 0) {
					opt->early_termination = 0;
				}
				else if (strcmp(long_options[option_index].name, "strike-box") == 0) {
					if (sscanf(optarg, "%lf", &(opt->strike_box)) != 1) {
						fprintf(stderr, "Invalid strike-box value\n");
						opt->strike_box = 1.0;
					}
				}
				else if (strcmp(long_options[option_index].name, "max-strikes") == 0) {
					if (sscanf(optarg, "%d", &(opt->max_strikes)) != 1) {
						fprintf(stderr, "Invalid max-strikes value\n");
						opt->max_strikes = 6;
					}
				}
				else if (strcmp(long_options[option_index].name, "enable-pruning") == 0) {
					opt->enable_pruning = 1;
				}
				else if (strcmp(long_options[option_index].name, "disable-pruning") == 0) {
					opt->enable_pruning = 0;
				}
				else if (strcmp(long_options[option_index].name, "pruning-factor") == 0) {
					if (sscanf(optarg, "%lf", &(opt->pruning_factor)) != 1) {
						fprintf(stderr, "Invalid pruning-factor value\n");
						opt->pruning_factor = 2.0;
					}
				}
				break;
			case 'h':
				print_help_statement();
				exit(0);
				break;
			case 'w':
				opt->use_nw=1;
				break;
			case 'p': //--paired
				strcpy(opt->paired_or_single,"paired");
				break;
			case 's': //--single
				strcpy(opt->paired_or_single,"single");
				break;
			case 'q': //--fastq
				opt->fastq=1;
				break;
			case '6':
				opt->skip_build=1;
				break;
			case '5': //print node info
				success = sscanf(optarg, "%s", opt->print_node_info);
				if (!success)
					fprintf(stderr, "Invalid file\n");
				break;
			case 'U':
				opt->unassigned=1;
				break;
			case 'r': //--reference
				opt->reference_mode = 1;
				break;
			case 'P': //--print-alignments
				opt->print_alignments = 1;
				break;
			case 'y':
				opt->use_partitions = 1;
				break;
			case 'v':
				opt->reverse_single_read=1;
				break;
			case '7':
				opt->print_all_nodes=1;
				break;
			case 'z':
				opt->reverse_second_of_paired_read=1;
				break;	
			case 'f':
				success = sscanf(optarg, "%s", opt->reference_file);
				if (!success)
					fprintf(stderr, "Invalid reference file.\n");
				break;
			case 'u':
				success = sscanf(optarg, "%lf", &(opt->score_constant));
				if (!success)
					fprintf(stderr, "Could not read score constant\n");
				break;
			case 't':
				success = sscanf(optarg, "%s", opt->print_trees_dir);
				if (!success)
					fprintf(stderr, "Invalid directory to print Newick trees.\n");
				break;
			case 'm':
				success = sscanf(optarg, "%s", opt->msa_file);
				if (!success)
					fprintf(stderr, "Invalid MSA file\n");
				break;
			case 'd':
				success = sscanf(optarg, "%s", opt->partitions_directory);
				if (!success)
					fprintf(stderr, "Invalid partitions output directory\n");
				break;
			case 'o':
				success = sscanf(optarg, "%s", opt->results_file);
				if (!success)
					fprintf(stderr, "Invalid output results file\n");
				break;
			case 'x':
				success = sscanf(optarg, "%s", opt->taxonomy_file);
				if (!success)
					fprintf(stderr, "Invalid taxonomy file\n");
				break;
			case 'g':
				success = sscanf(optarg, "%s", opt->read1_file);
				if (!success)
					fprintf(stderr, "Invalid single read file\n");
				break;
			case '1':
				success = sscanf(optarg, "%s", opt->read1_file);
				if (!success)
					fprintf(stderr, "Invalid read 1 file.\n");
				break;
			case '2':
				success = sscanf(optarg, "%s", opt->read2_file);
				if (!success)
					fprintf(stderr, "Invalid read 2 file.\n");
				break;
			case 'a':
				success = sscanf(optarg, "%s", opt->fasta_file);
				if (!success)
					fprintf(stderr, "Invalid fasta file.\n");
				break;
			case '3':
				success = sscanf(optarg, "%s", opt->print_alignments_dir);
				opt->print_alignments_to_file = 1;
				if (!success)
					fprintf(stderr, "Invalid directory.\n");
				break;
			case '4':
				success = sscanf(optarg, "%s", opt->treedir);
				if (!success)
					fprintf(stderr, "Invalid directory.\n");
				break;
			case 'c':
				success = sscanf(optarg, "%lf", &(opt->cinterval));
				if (!success)
					fprintf(stderr, "Invalid interval\n");
				break;
			case 'e':
				opt->use_leaf_portion = 1;
				break;
			case 'n':
				success = sscanf(optarg, "%d", &(opt->padding));
				if (!success)
					fprintf(stderr, "Invalid directory");
				break;
			case 'C':
				success = sscanf(optarg, "%d", &(opt->number_of_cores));
				if (!success)
					fprintf(stderr, "Invalid number");
				break;
			case 'L':
				success = sscanf(optarg, "%d", &(opt->number_of_lines_to_read));
				if (!success)
					fprintf(stderr, "Invalid number");
				break;
			case 'V':
				if (optarg) {
					success = sscanf(optarg, "%d", &(opt->verbose_level));
					if (!success || opt->verbose_level < 0 || opt->verbose_level > 3) {
						fprintf(stderr, "Invalid verbose level (0-3)\n");
						opt->verbose_level = 2;  // Default to INFO
					}
				} else {
					opt->verbose_level = 2;  // Default to INFO if no level specified
				}
				break;
			case 'l':
				strncpy(opt->log_file, optarg, sizeof(opt->log_file) - 1);
				opt->log_file[sizeof(opt->log_file) - 1] = '\0';
				break;
			case 'R':
				opt->enable_resource_monitoring = 1;
				break;
			case 'T':
				opt->enable_timing = 1;
				break;
			default:
				// Ignore unknown options, getopt will handle them
				break;
		}
	}
}
